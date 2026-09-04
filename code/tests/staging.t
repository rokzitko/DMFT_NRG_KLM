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

sub write_file {
    my ($path, $contents) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print $fh $contents;
    close($fh);
}

subtest "next-input alias migration" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("param.mu", "0.25\n");
    write_file("Delta.dat", "-1 0\n1 0\n");
    write_file("ReDelta.dat", "-1 0\n1 0\n");
    write_file("ImDelta.dat", "-1 -1\n1 -1\n");

    is(system($^X, "$scripts/next_aliases"), 0, "legacy inputs migrate");
    is(readlink("param.mu"), "param.mu.next", "param.mu is a next-input alias");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "ReDelta alias is explicit");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "ImDelta alias is explicit");
    is(readlink("Delta.dat"), "Delta.next.dat", "Gamma alias is explicit");
    ok(-s "param.mu.next", "mu target exists");
    is(system($^X, "$scripts/next_aliases"), 0, "alias maintenance is idempotent");
};

subtest "initialization treats Gamma as authoritative" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    install_numerical_stubs($dir);
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file(
        "param.loop",
        "clipDelta=0.1\nbroaden_max=2\nbroaden_ratio=2\nbroaden_min=0.5\n",
    );
    write_file("param.eps", "0.25\n");
    write_file("param.mu", "0\n");
    write_file("Delta.dat", "-2 0.1\n-1 0.2\n1 0.3\n2 0.1\n");
    write_file("ReDelta.dat", "-2 9\n-1 9\n1 9\n2 9\n");
    write_file("ImDelta.dat", "-2 -9\n-1 -9\n1 -9\n2 -9\n");

    is(system($^X, "$scripts/dmft_init"), 0,
       "preexisting Gamma restart initializes successfully");
    is(readlink("Delta.dat"), "Delta.next.dat", "Gamma alias is retained");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "real alias is regenerated");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "imaginary alias is regenerated");
    my @gamma = read_table("Delta.next.dat");
    my @real = read_table("ReDelta.next.dat");
    my @imaginary = read_table("ImDelta.next.dat");
    is($gamma[0][1], 0, "restart projection creates the lower support guard");
    is($gamma[-1][1], 0, "restart projection creates the upper support guard");
    my $consistent = 1;
    for my $index (0 .. $#gamma) {
        $consistent = 0 if abs($gamma[$index][1] + $imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "restart ImDelta is derived from Gamma");
    ok(!grep({ abs($_->[1] - 9) < 1e-15 } @real),
       "stale restart ReDelta values are discarded");
    is(scalar(@gamma), scalar(read_table("mesh.next.dat")),
       "restart mesh follows remeshed Gamma");
};

subtest "launcher migrates an Im-only legacy result" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("res", "scripts");
    install_numerical_stubs($dir);
    local $ENV{PATH} = "$dir/bin:$original_path";
    symlink("$scripts/migrate_legacy_res", "scripts/migrate_legacy_res") or die $!;
    symlink("$scripts/causalDelta", "scripts/causalDelta") or die $!;
    write_file("scripts/validate_bare_inputs", "#!/bin/sh\nexit 0\n");
    chmod(0755, "scripts/validate_bare_inputs");
    my $negative = "-2 -0.1\n-1 -0.2\n1 -0.2\n2 -0.1\n";
    my $positive = "-2 0\n-1 0.2\n1 0.2\n2 0\n";
    write_file("param.loop", "clipDelta=0.1\n");
    write_file("param.eps", "0.25\n");
    write_file("DOS.dat", "placeholder\n");
    write_file("CONVERGED", "\n");
    write_file("res/ImDelta.dat", $negative);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("occupancy.log", "0.1 0.15 0.05 0.7\n0.15 0.2 0.05 0.8\n");

    isnt(system($^X, "$scripts/DMFT"), 0,
         "launcher reaches the convergence guard after migrating legacy results");
    ok(!-e "res", "legacy result directory is removed");
    is(read_file("param.mu.used"), "0.15\n", "used mu is recovered from occupancy history");
    my @gamma = read_table("Delta.used.dat");
    my @imaginary = read_table("ImDelta.used.dat");
    is_deeply([map { $_->[1] } @gamma], [0, 0.2, 0.2, 0],
              "legacy ImDelta is canonicalized as authoritative Gamma");
    my $consistent = 1;
    for my $index (0 .. $#gamma) {
        $consistent = 0 if abs($gamma[$index][1] + $imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "legacy ImDelta is regenerated from Gamma");
    ok(-s "ReDelta.used.dat", "causal ReDelta is generated without a legacy real part");
    is(read_file("mesh.used.dat"), read_file("Delta.used.dat"),
       "used mesh follows authoritative Gamma");
    ok(-s "c-imG.dat" && -s "self.dat", "legacy result files move to root");
};

subtest "interrupted remeshing publication" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    install_numerical_stubs($dir);
    local $ENV{PATH} = "$dir/bin:$original_path";
    my $gamma = "-2 0\n-1 0.3\n1 0.3\n2 0\n";
    my $re = "-2 9\n-1 9\n1 9\n2 9\n";
    my $im = "-2 -8\n-1 -8\n1 -8\n2 -8\n";
    my $mesh = "-2 0\n-1 0\n1 0\n2 0\n";
    write_file("param.loop", "clipDelta=0.1\n");
    write_file("param.eps", "0\n");
    write_file("param.mu.next", "0\n");
    write_file("Delta.next.dat", "");
    write_file("ReDelta.next.dat", "old\n");
    write_file("ImDelta.next.dat", "old\n");
    write_file("Delta.next.dat.resampled", $gamma);
    write_file("ReDelta.next.dat.resampled", $re);
    write_file("ImDelta.next.dat.resampled", $im);
    write_file("mesh.next.dat.new", $mesh);
    write_file("mesh.next.dat.resample-ready", "ready\n");

    is(system($^X, "$scripts/dmft_init"), 0,
       "initialization recovers staged remeshing despite truncated live Gamma");
    my @restored_gamma = read_table("Delta.next.dat");
    my @restored_real = read_table("ReDelta.next.dat");
    my @restored_imaginary = read_table("ImDelta.next.dat");
    is_deeply([map { $_->[1] } @restored_gamma], [0, 0.3, 0.3, 0],
              "Gamma is restored and projected");
    ok(!grep({ abs($_->[1] - 9) < 1e-15 } @restored_real),
       "cached ReDelta is reconstructed instead of trusted");
    my $consistent = 1;
    for my $index (0 .. $#restored_gamma) {
        $consistent = 0
            if abs($restored_gamma[$index][1] + $restored_imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "cached ImDelta is reconstructed from Gamma");
    is(read_file("mesh.next.dat"), $mesh, "mesh is restored");
};

subtest "Broyden mixes Gamma with optional mu" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "res", "bin");
    my $mesh = "-2 0\n-1 0.2\n1 0.2\n2 0\n";
    my $raw = "-2 0\n-1 0.4\n1 0.4\n2 0\n";
    write_file("dmft/1-Delta.dat", $mesh);
    write_file("dmft/1-param.mu", "0\n");
    write_file("dmft/1-mu-occup.dat", "0 0.7\n");
    write_file("res/Delta.dat.OLD-common", $mesh);
    write_file("res/Delta.dat.NEW-common", $raw);
    write_file("res/param.mu.used", "0\n");
    write_file("param.eps", "0.25\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "alpha) printf 0.5;;\n" .
        "under) printf 0.5;;\n" .
        "goal) printf 0.8;;\n" .
        "maxdx) printf 0.05;;\n" .
        "broydenM) printf 50;;\n" .
        "clipDelta) printf 0.1;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    write_kk_stub("bin/kk");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--delta", "--mu",
        "--iteration", "1", "--workdir", "res",
        "--mu-input", "res/param.mu.used",
        "--mu-output", "res/param.mu.next",
        "--history-dir", "dmft", "--log", "res/occupancy.log",
    ), 0, "combined Gamma and mu update succeeds");
    my @mixed = read_table("res/Delta.dat.TEMP");
    cmp_ok(abs($mixed[1][1] - 0.3), "<=", 1e-15,
           "Broyden state contains one mixed Gamma component");
    cmp_ok(0 + read_file("res/param.mu.next"), ">", 0,
           "mu follows the Gamma block in the combined vector");

    is(system(
        $^X, "$scripts/causalDelta", "res/Delta.dat.TEMP",
        "res/Delta.next.dat", "res/ReDelta.next.dat", "res/ImDelta.next.dat",
    ), 0, "Broyden output passes through causal projection");
    my @projected = read_table("res/Delta.next.dat");
    is($projected[0][1], 0, "projection clears the lower raw endpoint");
    is($projected[-1][1], 0, "projection clears the upper raw endpoint");
    cmp_ok(abs($projected[1][1] - 0.3), "<=", 1e-15,
            "projection retains the mixed interior Gamma");
};

subtest "Broyden falls back when acceleration violates Gamma positivity" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "res", "bin");
    write_file("dmft/1-Delta.dat", "-2 0\n-1 0.1\n1 0.1\n2 0\n");
    write_file("dmft/1-Delta.raw.dat", "-2 0\n-1 0.111\n1 0.111\n2 0\n");
    write_file("dmft/2-Delta.dat", "-2 0\n-1 0.2\n1 0.2\n2 0\n");
    write_file("res/Delta.dat.OLD-common", "-2 0\n-1 0.2\n1 0.2\n2 0\n");
    write_file("res/Delta.dat.NEW-common", "-2 0\n-1 0.22\n1 0.22\n2 0\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "alpha) printf 0.5;;\n" .
        "broydenM) printf 50;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--delta", "--iteration", "2",
        "--workdir", "res", "--history-dir", "dmft",
    ), 0, "unstable accelerated proposal falls back successfully");
    my @mixed = read_table("res/Delta.dat.TEMP");
    cmp_ok(abs($mixed[1][1] - 0.21), "<=", 2e-15,
           "fallback uses the current linear Gamma step");
    cmp_ok(abs($mixed[2][1] - 0.21), "<=", 2e-15,
           "fallback applies consistently across Gamma components");
};

subtest "Broyden migrates legacy imaginary history" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "res", "bin");
    write_file("dmft/1-Delta.dat", "-2 0\n-1 0.2\n1 0.2\n2 0\n");
    write_file("dmft/1-ImDelta.raw.dat", "-2 -0.1\n-1 -0.4\n1 -0.4\n2 -0.1\n");
    write_file("dmft/2-Delta.dat", "-2 0\n-1 0.25\n1 0.25\n2 0\n");
    write_file("res/Delta.dat.OLD-common", "-2 0\n-1 0.25\n1 0.25\n2 0\n");
    write_file("res/Delta.dat.NEW-common", "-2 0\n-1 0.5\n1 0.5\n2 0\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "alpha) printf 0.5;;\n" .
        "broydenM) printf 50;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--delta", "--iteration", "2",
        "--workdir", "res", "--history-dir", "dmft",
    ), 0, "Gamma Broyden accepts converted legacy history");
    my @migrated = read_table("dmft/1-Delta.raw.dat");
    is_deeply([map { $_->[1] } @migrated], [0, 0.4, 0.4, 0],
              "legacy raw ImDelta is migrated with canonical endpoint guards");
    ok(-s "res/Delta.dat.TEMP", "migrated history contributes to a Gamma update");
    my @mixed = read_table("res/Delta.dat.TEMP");
    is($mixed[0][1], 0, "Broyden excludes the lower endpoint from its state");
    is($mixed[-1][1], 0, "Broyden excludes the upper endpoint from its state");
};

subtest "Broyden preserves compatibility mu symlink" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "bin");
    write_file("param.mu.next", "0\n");
    symlink("param.mu.next", "param.mu") or die $!;
    write_file("dmft/1-param.mu", "0\n");
    write_file("dmft/1-mu-occup.dat", "0 0.7\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "under) printf 0.5;;\n" .
        "goal) printf 0.8;;\n" .
        "maxdx) printf 0.05;;\n" .
        "broydenM) printf 50;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--mu", "--iteration", "1",
    ), 0, "Broyden mu update succeeds through the legacy name");
    is(readlink("param.mu"), "param.mu.next", "Broyden does not replace the symlink");
    cmp_ok(0 + read_file("param.mu.next"), ">", 0, "Broyden updates the explicit target");
};

subtest "accurate occupancy is staged before convergence and reuses H0" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("1", "bin", "scripts");
    for my $name (qw(
        bandDOS occupancy_control dmftDOS-stable causalDelta copyresults next_aliases
    )) {
        symlink("$scripts/$name", "scripts/$name") or die $!;
    }

    my $mesh = "-1 0\n-0.5 0.2\n0.5 0.2\n1 0\n";
    my $sigma_real = "-1 0\n-0.5 0\n0.5 0\n1 0\n";
    my $sigma_imaginary = "-1 -0.1\n-0.5 -0.1\n0.5 -0.1\n1 -0.1\n";
    write_file("param.loop", "placeholder=true\n");
    write_file("param.eps", "0\n");
    write_file("param.mu.next", "0\n");
    write_file("DOS.dat", $mesh);
    write_file("1/DONE", "\n");

    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "Nz) printf 1;; mixing) printf linear;; mucontrol) printf std;;\n" .
        "occupancy_mode) printf accurate;; track) printf false;;\n" .
        "clipDelta) printf 0.1;; clipSigma) printf 1e-12;;\n" .
        "T) printf 0.1;; goal) printf 0.8;; maxdx) printf 0.2;; under) printf 0.5;;\n" .
        "occupancy_solve_tol) printf 1e-8;; occupancy_mu_tol) printf 1e-10;;\n" .
        "occupancy_weight_tol) printf 1e-4;; occupancy_integ_epsabs) printf 1e-10;;\n" .
        "occupancy_integ_epsrel) printf 1e-9;; occupancy_maxeval) printf 12;;\n" .
        "alpha) printf 0.5;; opt_delta) printf 1e-3;; opt_mindx) printf 1e-10;;\n" .
        "opt_ymin) printf 1e-5;; *) exit 1;; esac\n",
    );
    write_file(
        "bin/hilb",
        "#!/usr/bin/env perl\nuse strict; use warnings;\n" .
        "my (\$moment, \$mu) = (0, 0);\n" .
        "for (my \$i = 0; \$i < \@ARGV; ++\$i) {\n" .
        "  \$moment = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-n';\n" .
        "  \$mu = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-x';\n" .
        "}\n" .
        "open(my \$events, '>>', \$ENV{EVENTS}) or die \$!;\n" .
        "print {\$events} \"hilb \$moment \$mu\\n\"; close(\$events);\n" .
        "my (\$real, \$imaginary) = \@ARGV[-2, -1];\n" .
        "open(my \$rf, '>', \$real) or die \$!;\n" .
        "open(my \$if, '>', \$imaginary) or die \$!;\n" .
        "if (\$moment == 1) {\n" .
        "  print {\$rf} \"-1 0\\n-0.5 0\\n0.5 0\\n1 0\\n\";\n" .
        "  print {\$if} \"-1 0\\n-0.5 -1\\n0.5 -1\\n1 0\\n\";\n" .
        "} else {\n" .
        "  my \$marker = 0.3 + \$mu;\n" .
        "  print {\$rf} \"-1 1\\n-0.5 1\\n0.5 1\\n1 1\\n\";\n" .
        "  print {\$if} \"-1 \$marker\\n-0.5 0.2\\n0.5 0.2\\n1 0.1\\n\";\n" .
        "}\nclose(\$rf); close(\$if);\n",
    );
    write_file(
        "bin/integ",
        "#!/usr/bin/env perl\nuse strict; use warnings;\n" .
        "my \$file = \$ARGV[-1]; open(my \$fh, '<', \$file) or die \$!;\n" .
        "my (\$x, \$marker) = split(/\\s+/, <\$fh>); close(\$fh);\n" .
        "if (grep { \$_ eq '--total' } \@ARGV) { print \"1\\n\"; }\n" .
        "else { print \$marker, \"\\n\"; }\n",
    );
    write_file(
        "scripts/average",
        "#!/bin/sh\nfor name in c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat; do\n" .
        "  printf '%s\\n' '-1 0.3' '-0.5 0.2' '0.5 0.2' '1 0.1' >\"\$1/\$name\"\n" .
        "done\n",
    );
    write_file("scripts/realparts", "#!/bin/sh\nexit 0\n");
    write_file(
        "scripts/sigmatrick",
        "#!/bin/sh\nprintf '%s' '$sigma_real' >\"\$1/resigma.dat\"\n" .
        "printf '%s' '$sigma_imaginary' >\"\$1/imsigma.dat\"\n" .
        "printf '%s\\n' '-1 0.3' '-0.5 0.2' '0.5 0.2' '1 0.1' >\"\$1/c-self.dat\"\n",
    );
    write_file(
        "scripts/diffs",
        "#!/bin/sh\nprintf '%s\\n' '1 1e-3' >res/DIFFS_C\n" .
        "printf '%s\\n' '1 1e-3' >res/DIFFS_LatLoc\n",
    );
    write_file(
        "scripts/checkconv",
        "#!/bin/sh\ngrep -q '^mode=accurate\$' res/OCCUPANCY_METRICS || exit 81\n" .
        "test -s res/H0.re.dat -a -s res/H0.im.dat -a -s res/H0.mu || exit 82\n" .
        "printf '%s\\n' checkconv >>\"\$EVENTS\"\n" .
        "while [ \"\$#\" -gt 0 ]; do\n" .
        "  if [ \"\$1\" = --decision ]; then printf '%s\\n' stop >\"\$2\"; exit 0; fi\n" .
        "  shift\ndone\nexit 83\n",
    );
    write_file("scripts/DMFT", "#!/bin/sh\nexit 91\n");
    write_file("bin/gatherlastlines", "#!/bin/sh\nprintf '%s\\n' value >\"\$1\"\n");
    write_file("bin/columnavg_comment", "#!/bin/sh\nprintf '%s\\n' value\n");
    write_file(
        "bin/optimize_mesh",
        "#!/usr/bin/env perl\nopen(my \$fh, '<', \$ARGV[0]) or die \$!; print while <\$fh>;\n",
    );
    write_file(
        "bin/mesh_union",
        "#!/usr/bin/env perl\nopen(my \$fh, '<', \$ARGV[0]) or die \$!; print while <\$fh>;\n",
    );
    write_file(
        "bin/resample",
        "#!/usr/bin/env perl\nuse File::Copy qw(copy); copy(\$ARGV[-3], \$ARGV[-1]) or die \$!;\n",
    );
    write_file(
        "bin/mixy",
        "#!/usr/bin/env perl\nopen(my \$fh, '<', \$ARGV[1]) or die \$!; print while <\$fh>;\n",
    );
    write_file("bin/delete_NRG_output", "#!/bin/sh\nexit 0\n");
    write_kk_stub("bin/kk");
    chmod(0755, glob("bin/*"), glob("scripts/*"));

    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";
    local $ENV{EVENTS} = "$dir/events";
    is(system(
        $^X, "scripts/causalDelta", "DOS.dat",
        "Delta.next.dat", "ReDelta.next.dat", "ImDelta.next.dat",
    ), 0, "accurate-cycle fixture starts from a causal Gamma triplet");
    write_file("mesh.next.dat", $mesh);

    is(system($^X, "$scripts/dmft_done"), 0,
       "one complete accurate standard-control cycle stages and publishes");
    ok(!-e "res", "successful cycle removes its staging directory");
    ok(-e "STOP", "stub convergence decision is published");
    cmp_ok(abs((0.0 + read_file("param.mu.next")) - 0.05), "<=", 2e-12,
           "published next mu is the underrelaxed occupancy update");
    like(read_file("events"), qr/checkconv\nhilb 1 0\.05[^\n]*\n\z/,
         "convergence sees occupancy/cache before dmftDOS evaluates only H1");
    like(read_file("OCCUPANCY_METRICS"), qr/^mode=accurate$/m,
         "accurate occupancy metrics publish with the result batch");
};

subtest "ready publication and cleanup" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("res", "scripts", "bin");
    symlink("$scripts/copyresults", "scripts/copyresults") or die $!;
    symlink("$scripts/next_aliases", "scripts/next_aliases") or die $!;

    write_file("scripts/DMFT", "#!/bin/sh\nexit 0\n");
    write_file(
        "scripts/post",
        "#!/bin/sh\n" .
        "[ \"\$(cat param.mu.used)\" = 0.1 ] || exit 91\n" .
        "[ \"\$(cat DIFFS_C)\" = '1 1e-12' ] || exit 92\n" .
        "grep -q '^n_old=0.8\$' OCCUPANCY_METRICS || exit 93\n" .
        "printf '%s\\n' published >post.seen\n",
    );
    write_file("bin/delete_NRG_output", "#!/bin/sh\nexit 0\n");
    write_getparam_stub("bin/getparam");
    write_kk_stub("bin/kk");
    chmod(0755, "scripts/DMFT", "scripts/post", glob("bin/*"));
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    my $negative = "-2 0\n-1 -0.2\n1 -0.2\n2 0\n";
    my $positive = "-2 0\n-1 0.2\n1 0.2\n2 0\n";
    my $mesh = "-2 0\n-1 0\n1 0\n2 0\n";
    write_file("param.loop", "track=false\nclipDelta=0.1\n");
    write_file("param.eps", "0\n");
    write_causal_stage_triplet("used", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_causal_stage_triplet("next", $positive);
    write_file("res/param.mu.next", "0.2\n");
    write_file("res/mesh.next.dat", $mesh);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("res/$_", "value\n") for qw(custom custom.avg customfdm customfdm.avg);
    write_file("res/DIFFS_C", "1 1e-12\n");
    write_file("res/occupancy.log", "0.1 0.2 0.1 0.8\n");
    write_file("res/OCCUPANCY_METRICS", "version=1\niteration=1\nn_old=0.8\n");
    write_file("res/ITER", "1\n");
    write_file("res/DECISION", "converged\n");
    write_file("res/READY", "source_iter=0\ntarget_iter=1\ndecision=converged\n");

    write_file("res/ReDelta.next.dat", "-2 9\n-1 9\n1 9\n2 9\n");
    isnt(system($^X, "$scripts/dmft_done", "--publish"), 0,
         "ready replay rejects a stale independent ReDelta");
    ok(-d "res", "rejected ready transaction remains available for repair");
    write_causal_stage_triplet("next", $positive);

    is(system($^X, "$scripts/dmft_done", "--publish"), 0, "ready stage publishes");
    ok(!-e "res", "staging directory is removed");
    is(readlink("param.mu"), "param.mu.next", "compatibility mu alias is installed");
    is(readlink("mesh.dat"), "mesh.next.dat", "compatibility mesh alias is installed");
    is(read_file("ITER"), "1\n", "iteration is published");
    is(read_file("param.mu.used"), "0.1\n", "used mu remains with results");
    is(read_file("param.mu.next"), "0.2\n", "next mu remains distinct");
    ok(-s "c-imG.dat" && -s "self.dat", "result batch is in the main directory");
    is(read_file("post.seen"), "published\n",
       "terminal postprocessing sees current diagnostics and used mu");
    ok(!-e "POSTPROCESSING_FAILED", "diagnostic publication completed before post");
    is(read_file("OCCUPANCY_METRICS"), "version=1\niteration=1\nn_old=0.8\n",
       "occupancy metrics are published with the result batch");

    make_path("res");
    write_causal_stage_triplet("used", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_causal_stage_triplet("next", $positive);
    write_file("res/param.mu.next", "0.2\n");
    write_file("res/mesh.next.dat", $mesh);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("res/$_", "value\n") for qw(custom custom.avg customfdm customfdm.avg);
    write_file("res/ITER", "1\n");
    write_file("res/DECISION", "continue\n");
    write_file("res/READY.complete", "source_iter=0\ntarget_iter=1\ndecision=continue\n");
    write_file("ITER", "x");

    is(system($^X, "$scripts/dmft_done", "--publish"), 0,
        "backup marker repairs malformed ITER during replay");
    is(read_file("ITER"), "1\n", "replay restores the iteration file");
    ok(!-e "res", "replayed staging directory is removed");

    make_path("res");
    write_file("res/READY.complete", "source_iter=0\ntarget_iter=1\ndecision=continue\n");
    write_file("res/partially-removed-file", "leftover\n");
    is(system($^X, "$scripts/dmft_done", "--publish"), 0,
        "published iteration finishes interrupted staging cleanup");
    ok(!-e "res", "partially removed staging directory is cleaned");
};

chdir($original) or die $!;
done_testing();

sub read_file {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    local $/;
    my $contents = <$fh>;
    close($fh);
    return $contents;
}

sub read_table {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    my @rows;
    while (<$fh>) {
        next if /^\s*$/;
        my @fields = split;
        @fields == 2 or die "Expected two columns in $path at line $.";
        push(@rows, [map { 0.0 + $_ } @fields]);
    }
    close($fh) or die "Can't close $path: $!";
    return @rows;
}

sub write_causal_stage_triplet {
    my ($kind, $gamma) = @_;
    my $raw = "res/.Delta.$kind.raw.$$";
    write_file($raw, $gamma);
    system(
        $^X, "$scripts/causalDelta", $raw,
        "res/Delta.$kind.dat", "res/ReDelta.$kind.dat", "res/ImDelta.$kind.dat",
    ) == 0 or die "Failed to create causal $kind staging fixture: $?";
    unlink($raw) or die "Can't remove $raw: $!";
}

sub install_numerical_stubs {
    my ($directory) = @_;
    make_path("$directory/bin");
    write_getparam_stub("$directory/bin/getparam");
    write_kk_stub("$directory/bin/kk");
    write_file("$directory/bin/mesh_union", <<'MESH_UNION');
#!/usr/bin/env perl
use strict;
use warnings;
my $source = $ARGV[-1];
open(my $fh, "<", $source) or die $!;
print while <$fh>;
close($fh) or die $!;
MESH_UNION
    write_file("$directory/bin/resample", <<'RESAMPLE');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy qw(copy);
copy($ARGV[-3], $ARGV[-1]) or die $!;
RESAMPLE
    chmod(0755, glob("$directory/bin/*"));
}

sub write_getparam_stub {
    my ($path) = @_;
    write_file($path, <<'GETPARAM');
#!/usr/bin/env perl
use strict;
use warnings;
my ($name, $file) = @ARGV;
defined($name) && defined($file) or exit 1;
open(my $fh, "<", $file) or exit 1;
while (<$fh>) {
    if (/^\s*\Q$name\E\s*=\s*(.*?)\s*$/) {
        print "$1\n";
        exit 0;
    }
}
exit 1;
GETPARAM
    chmod(0755, $path);
}

sub write_kk_stub {
    my ($path) = @_;
    write_file($path, <<'KK');
#!/usr/bin/env perl
use strict;
use warnings;
my ($input, $output) = @ARGV[-2, -1];
open(my $in, "<", $input) or die $!;
open(my $out, ">", $output) or die $!;
while (<$in>) {
    next if /^\s*(?:#.*)?$/;
    my ($omega) = split;
    defined($omega) or die "invalid input\n";
    print {$out} "$omega 0\n";
}
close($in) or die $!;
close($out) or die $!;
KK
    chmod(0755, $path);
}
