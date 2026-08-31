#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use Test::More;

my $scripts = "$FindBin::Bin/../scripts";
my $original = getcwd();
my $original_path = $ENV{PATH};

subtest "fast occupancy control shifts frequencies without extrapolation" => sub {
    my $dir = occupancy_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{INTEG_LOG} = "$dir/integ.log";

    ok(-x "$scripts/occupancy_control", "occupancy controller is executable");
    is(
        system(
            "$scripts/occupancy_control", "--fast",
            "--spectrum", "imaw.dat",
            "--mu-input", "mu.used",
            "--mu-output", "mu.next",
            "--log", "occupancy.log",
            "--metrics", "metrics",
            "--iteration", "3",
        ),
        0,
        "fast controller finds a safeguarded update",
    );

    cmp_ok(abs(read_scalar("mu.next") - 0.05), "<=", 2e-12,
           "underrelaxation is applied to the bracketed root");
    my @log = split(/\s+/, read_file("occupancy.log"));
    is(scalar(@log), 4, "legacy occupancy log retains four columns");
    cmp_ok(abs($log[0]), "<=", 1e-15, "log records the used mu");
    cmp_ok(abs($log[1] - 0.05), "<=", 2e-12, "log records the next mu");
    cmp_ok(abs($log[2] - 0.05), "<=", 2e-12, "log records the applied step");
    cmp_ok(abs($log[3] - 0.6), "<=", 2e-12, "log records current occupancy");

    my %metrics = read_metrics("metrics");
    is($metrics{mode}, "fast", "metrics identify the fast evaluator");
    is($metrics{status}, "bracketed", "metrics identify a bracketed solve");
    cmp_ok(abs($metrics{mu_root} - 0.1), "<=", 2e-12,
           "the undamped filling root is recorded");
    cmp_ok(abs($metrics{n_new} - 0.7), "<=", 2e-12,
           "final occupancy is evaluated at the damped mu");
    cmp_ok(abs($metrics{convergence_dx} - 0.05), "<=", 2e-12,
           "proposed mu movement is available to convergence checks");

    my @integrations = split(/\n/, read_file("integ.log"));
    my @first_x = map { /first_x=([^\s]+)/ ? 0.0 + $1 : () } @integrations;
    ok(grep(abs($_ + 1.4) <= 2e-12, @first_x),
       "bracketing trial shifts the mesh left when mu increases");
    ok(grep(abs($_ + 1.05) <= 2e-12, @first_x),
       "damped final trial shifts only the frequency column");
};

subtest "accurate occupancy control evaluates frozen Sigma and publishes H0" => sub {
    my $dir = occupancy_fixture();
    chdir($dir) or die $!;
    make_path("res", "dmft");
    write_file("res/imaw.dat", spectrum(0.3));
    write_file("res/reaw.dat", real_spectrum());
    write_file("res/resigma.dat", "-1 0\n0 0\n1 0\n");
    write_file("res/imsigma.dat", "-1 -0.1\n0 -0.1\n1 -0.1\n");
    write_file("DOS.dat", "-1 0\n0 1\n1 0\n");
    write_banddos_stub("bin/bandDOS");

    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{INTEG_LOG} = "$dir/integ.log";
    local $ENV{BANDDOS_LOG} = "$dir/banddos.log";
    is(
        system(
            "$scripts/occupancy_control", "--accurate",
            "--spectrum", "res/imaw.dat",
            "--real-spectrum", "res/reaw.dat",
            "--result-dir", "res",
            "--dos", "DOS.dat",
            "--mu-input", "mu.used",
            "--mu-output", "mu.next",
            "--log", "occupancy.log",
            "--metrics", "metrics",
            "--iteration", "4",
            "--banddos", "bin/bandDOS",
            "--cache-real", "res/H0.re.dat",
            "--cache-imaginary", "res/H0.im.dat",
            "--cache-mu", "res/H0.mu",
        ),
        0,
        "accurate controller solves with frozen-Sigma lattice evaluations",
    );

    my %metrics = read_metrics("metrics");
    is($metrics{mode}, "accurate", "metrics identify the accurate evaluator");
    my %cache_metadata = read_metrics("res/H0.mu");
    cmp_ok(abs($cache_metadata{mu} - read_scalar("mu.next")), "<=", 1e-15,
           "cache marker matches the published next mu");
    like($cache_metadata{resigma_sha256}, qr/\A[0-9a-f]{64}\z/,
         "cache marker identifies the frozen self-energy");
    cmp_ok(abs((read_table("res/H0.im.dat"))[0][1] - 0.35), "<=", 2e-12,
           "cached H0 spectrum is the final damped-mu evaluation");
    is(scalar(read_table("res/H0.re.dat")), 5,
       "matching real H0 component is published");

    my @trial_mu = map { 0.0 + $_ } grep { length } split(/\n/, read_file("banddos.log"));
    ok(grep(abs($_ - 0.4) <= 2e-12, @trial_mu),
       "accurate evaluator computes the bracket edge");
    ok(grep(abs($_ - 0.1) <= 2e-12, @trial_mu),
       "accurate evaluator computes the occupancy root");
    ok(grep(abs($_ - 0.05) <= 2e-12, @trial_mu),
       "accurate evaluator computes the underrelaxed final H0");

    write_file("occupancy.log", "-0.0125 0 0.0125 0.59\n");
    write_file("mu.unchanged", "7\n");
    unlink("banddos.log") or die $!;
    local $ENV{BANDDOS_FAIL} = 1;
    is(
        system(
            "$scripts/occupancy_control", "--accurate", "--measure-only", "--history",
            "--spectrum", "res/imaw.dat",
            "--mu-input", "mu.used",
            "--mu-output", "mu.unchanged",
            "--log", "occupancy.log",
            "--metrics", "measured.metrics",
            "--iteration", "5",
            "--history-dir", "dmft",
            "--banddos", "bin/bandDOS",
        ),
        0,
        "accurate measure-only mode does not invoke the Hilbert backend",
    );
    ok(!-e "banddos.log", "measure-only mode made no trial lattice evaluation");
    is(read_file("mu.unchanged"), "7\n", "measure-only mode preserves mu output");
    is(read_file("dmft/5-mu-occup.dat"), "0 0.59999999999999998\n",
       "Broyden history records the used mu and measured occupancy");
    my %measured = read_metrics("measured.metrics");
    is($measured{status}, "measured", "measurement status is explicit");
    cmp_ok(abs($measured{convergence_dx} - 0.0125), "<=", 1e-15,
           "measurement carries the previously applied Broyden step");

    write_file("occupancy.log", "0 0.5 0.5 0.59\n");
    write_file("stale.metrics", "preserve\n");
    isnt(
        system(
            "$scripts/occupancy_control", "--accurate", "--measure-only",
            "--spectrum", "res/imaw.dat", "--mu-input", "mu.used",
            "--mu-output", "mu.unchanged", "--log", "occupancy.log",
            "--metrics", "stale.metrics", "--iteration", "6",
            "--banddos", "bin/bandDOS",
        ),
        0,
        "measure-only mode rejects a log step not leading to the current mu",
    );
    is(read_file("stale.metrics"), "preserve\n",
       "stale Broyden log does not replace convergence metrics");

    write_file("occupancy.log", "0 0 0 0.6\n");
    write_file("failed.metrics", "preserve metrics\n");
    my $cache_before_failure = read_file("res/H0.mu");
    isnt(
        system(
            "$scripts/occupancy_control", "--accurate",
            "--spectrum", "res/imaw.dat", "--real-spectrum", "res/reaw.dat",
            "--result-dir", "res", "--dos", "DOS.dat",
            "--mu-input", "mu.used", "--mu-output", "mu.unchanged",
            "--log", "occupancy.log", "--metrics", "failed.metrics",
            "--iteration", "6", "--banddos", "bin/bandDOS",
            "--cache-real", "res/H0.re.dat",
            "--cache-imaginary", "res/H0.im.dat", "--cache-mu", "res/H0.mu",
        ),
        0,
        "trial Hilbert failure aborts an accurate update",
    );
    is(read_file("mu.unchanged"), "7\n", "trial failure preserves mu output");
    is(read_file("failed.metrics"), "preserve metrics\n",
       "trial failure preserves occupancy metrics");
    is(read_file("res/H0.mu"), $cache_before_failure,
       "trial failure preserves the previous complete H0 cache marker");
};

subtest "controller limits steps and preserves outputs on validation failure" => sub {
    my $dir = occupancy_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{INTEG_LOG} = "$dir/integ.log";
    local $ENV{GOAL} = 1.8;
    local $ENV{MAXDX} = 0.05;

    is(
        system(
            "$scripts/occupancy_control", "--fast",
            "--spectrum", "imaw.dat", "--mu-input", "mu.used",
            "--mu-output", "mu.next", "--log", "occupancy.log",
            "--metrics", "metrics", "--iteration", "2",
        ),
        0,
        "unbracketed target takes a bounded step",
    );
    cmp_ok(abs(read_scalar("mu.next") - 0.05), "<=", 1e-15,
           "limited update respects maxdx");
    my %limited = read_metrics("metrics");
    is($limited{status}, "limited", "limited status is recorded");

    write_file("mu.next", "9\n");
    write_file("occupancy.log", "old log\n");
    write_file("metrics", "old metrics\n");
    local $ENV{BAD_WEIGHT} = 1;
    isnt(
        system(
            "$scripts/occupancy_control", "--fast",
            "--spectrum", "imaw.dat", "--mu-input", "mu.used",
            "--mu-output", "mu.next", "--log", "occupancy.log",
            "--metrics", "metrics", "--iteration", "2",
        ),
        0,
        "spectral sum-rule failure is fatal",
    );
    is(read_file("mu.next"), "9\n", "failed update preserves chemical potential");
    is(read_file("occupancy.log"), "old log\n", "failed update preserves legacy log");
    is(read_file("metrics"), "old metrics\n", "failed update preserves metrics");
};

subtest "convergence requires filling and chemical-potential closure" => sub {
    my $dir = occupancy_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file("DIFFS_C", "11 1e-12\n");

    write_metrics_file("metrics", 11, 0.1, 0.0);
    is(run_checkconv(), 0, "large filling error is a valid nonconverged decision");
    is(read_file("DECISION"), "continue\n", "filling error blocks convergence");

    write_metrics_file("metrics", 11, 1e-10, 1e-3);
    is(run_checkconv(), 0, "large mu step is a valid nonconverged decision");
    is(read_file("DECISION"), "continue\n", "mu movement blocks convergence");

    write_metrics_file("metrics", 11, 1e-10, 1e-12);
    is(run_checkconv(), 0, "closed occupancy produces a decision");
    is(read_file("DECISION"), "converged\n",
       "spectral, filling, and mu criteria jointly converge");

    write_metrics_file("metrics", 10, 0.0, 0.0);
    write_file("DECISION", "preserve\n");
    isnt(run_checkconv(), 0, "stale occupancy metrics are rejected");
    is(read_file("DECISION"), "preserve\n", "invalid metrics preserve prior decision");
};

chdir($original) or die $!;
done_testing();

sub occupancy_fixture {
    my $dir = tempdir(CLEANUP => 1);
    make_path("$dir/bin");
    write_file("$dir/param.loop", "placeholder=true\n");
    write_file("$dir/mu.used", "0\n");
    write_file("$dir/imaw.dat", spectrum(0.3));
    write_file(
        "$dir/bin/getparam",
        "#!/bin/sh\n" .
        "case \"\$1\" in\n" .
        "T) printf '%s\\n' 0.1;;\n" .
        "goal) printf '%s\\n' \"\${GOAL:-0.8}\";;\n" .
        "maxdx) printf '%s\\n' \"\${MAXDX:-0.2}\";;\n" .
        "under) printf '%s\\n' 0.5;;\n" .
        "occupancy_mode) printf '%s\\n' accurate;;\n" .
        "occupancy_solve_tol) printf '%s\\n' 1e-8;;\n" .
        "occupancy_mu_tol) printf '%s\\n' 1e-10;;\n" .
        "occupancy_weight_tol) printf '%s\\n' 1e-4;;\n" .
        "occupancy_integ_epsabs) printf '%s\\n' 1e-10;;\n" .
        "occupancy_integ_epsrel) printf '%s\\n' 1e-9;;\n" .
        "occupancy_maxeval) printf '%s\\n' 12;;\n" .
        "clipSigma) printf '%s\\n' 1e-12;;\n" .
        "conveps) printf '%s\\n' 1e-8;;\n" .
        "convwindow) printf '%s\\n' 1;;\n" .
        "miniter) printf '%s\\n' 1;;\n" .
        "maxiter) printf '%s\\n' 100;;\n" .
        "*) exit 1;;\n" .
        "esac\n",
    );
    write_file(
        "$dir/bin/integ",
        "#!/usr/bin/env perl\n" .
        "use strict; use warnings;\n" .
        "my \$file = \$ARGV[-1];\n" .
        "open(my \$fh, '<', \$file) or die \$!;\n" .
        "my \$line = <\$fh>; close(\$fh);\n" .
        "my (\$x, \$marker) = split(/\\s+/, \$line);\n" .
        "my \$occupancy = 2 * \$marker + 2 * (-1 - \$x);\n" .
        "if (defined \$ENV{INTEG_LOG}) {\n" .
        "  open(my \$log, '>>', \$ENV{INTEG_LOG}) or die \$!;\n" .
        "  print {\$log} join(' ', \@ARGV), \" first_x=\$x occupancy=\$occupancy\\n\";\n" .
        "  close(\$log);\n" .
        "}\n" .
        "if (grep { \$_ eq '--total' } \@ARGV) {\n" .
        "  print \$ENV{BAD_WEIGHT} ? \"0.8\\n\" : \"1\\n\";\n" .
        "} else { print \$occupancy / 2, \"\\n\"; }\n",
    );
    chmod(0755, "$dir/bin/getparam", "$dir/bin/integ");
    return $dir;
}

sub write_banddos_stub {
    my ($path) = @_;
    write_file(
        $path,
        "#!/usr/bin/env perl\n" .
        "use strict; use warnings;\n" .
        "exit 81 if \$ENV{BANDDOS_FAIL};\n" .
        "my (\$mu, \$real, \$imaginary);\n" .
        "for (my \$i = 0; \$i < \@ARGV; ++\$i) {\n" .
        "  \$mu = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--mu-value';\n" .
        "  \$real = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--re-output';\n" .
        "  \$imaginary = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--im-output';\n" .
        "}\n" .
        "open(my \$log, '>>', \$ENV{BANDDOS_LOG}) or die \$!;\n" .
        "print {\$log} \"\$mu\\n\"; close(\$log);\n" .
        "open(my \$rf, '>', \$real) or die \$!;\n" .
        "print {\$rf} \"-1 -0.2\\n-0.5 -0.1\\n0 0\\n0.5 0.1\\n1 0.2\\n\"; close(\$rf);\n" .
        "open(my \$if, '>', \$imaginary) or die \$!;\n" .
        "my \$marker = 0.3 + \$mu;\n" .
        "print {\$if} \"-1 \$marker\\n-0.5 0.2\\n0 0.2\\n0.5 0.2\\n1 0.1\\n\"; close(\$if);\n",
    );
    chmod(0755, $path);
}

sub spectrum {
    my ($marker) = @_;
    return "-1 $marker\n-0.5 0.2\n0 0.2\n0.5 0.2\n1 0.1\n";
}

sub real_spectrum {
    return "-1 -0.2\n-0.5 -0.1\n0 0\n0.5 0.1\n1 0.2\n";
}

sub run_checkconv {
    return system(
        $^X, "$scripts/checkconv", "--iteration", "11",
        "--diffs", "DIFFS_C", "--occupancy-metrics", "metrics",
        "--decision", "DECISION",
    );
}

sub write_metrics_file {
    my ($path, $iteration, $error, $step) = @_;
    write_file(
        $path,
        "version=1\niteration=$iteration\nerror_old=$error\nconvergence_dx=$step\n",
    );
}

sub read_metrics {
    my ($path) = @_;
    my %values;
    for my $line (split(/\n/, read_file($path))) {
        next if $line eq "";
        my ($key, $value) = split(/=/, $line, 2);
        $values{$key} = $value;
    }
    return %values;
}

sub read_scalar {
    my ($path) = @_;
    my $value = read_file($path);
    $value =~ s/^\s+|\s+$//g;
    return 0.0 + $value;
}

sub write_file {
    my ($path, $contents) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print {$fh} $contents;
    close($fh) or die "Can't close $path: $!";
}

sub read_file {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    local $/;
    my $contents = <$fh>;
    close($fh) or die "Can't close $path: $!";
    return $contents;
}

sub read_table {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    my @rows;
    while (<$fh>) {
        next if /^\s*$/;
        my @fields = split;
        push(@rows, [map { 0.0 + $_ } @fields]);
    }
    close($fh) or die "Can't close $path: $!";
    return @rows;
}
