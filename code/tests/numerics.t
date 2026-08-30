#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use Test::More;

my $code = "$FindBin::Bin/..";
my $scripts = "$code/scripts";
my $original = getcwd();
my $original_path = $ENV{PATH};

require "$scripts/Bubble.pm";

subtest "Bethe tables use one normalized edge-refined mesh" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;

    local $ENV{PATH} = $original_path;
    my $old_umask = umask(0022);
    is(system("$code/mkDOS"), 0, "mkDOS accepts the analytic Steffen backend");
    is(system("$code/mkPHI"), 0, "mkPHI derives Phi from the normalized DOS");
    umask($old_umask);

    my @dos = read_table("DOS.dat");
    my @phi = read_table("PHI.dat");
    is(scalar(@dos), 2601, "DOS has the expected number of rows");
    is(scalar(@phi), 2601, "Phi has the expected number of rows");
    is(read_file("DOS.dat"), read_file("$code/DOS.dat"),
       "checked-in DOS matches the generator");
    is(read_file("PHI.dat"), read_file("$code/PHI.dat"),
       "checked-in Phi matches the generator");
    is((stat("DOS.dat"))[2] & 0777, 0644, "DOS uses the umask-derived file mode");
    is((stat("PHI.dat"))[2] & 0777, 0644, "Phi uses the umask-derived file mode");

    is($dos[300][0], -1, "lower band edge is exact");
    is($dos[2300][0], 1, "upper band edge is exact");
    is($dos[300][1], 0, "DOS vanishes at the lower edge");
    is($dos[2300][1], 0, "DOS vanishes at the upper edge");
    cmp_ok($dos[301][0] - $dos[300][0], "<",
           $dos[1301][0] - $dos[1300][0], "mesh is refined at band edges");

    my $same_mesh = 1;
    my $maximum_relation_error = 0;
    for my $index (0 .. $#dos) {
        $same_mesh = 0 if $phi[$index][0] != $dos[$index][0];
        my $expected = $dos[$index][1]
            * max(0, 1 - $dos[$index][0] ** 2);
        $maximum_relation_error = max(
            $maximum_relation_error,
            abs($phi[$index][1] - $expected),
        );
    }
    ok($same_mesh, "DOS and Phi use exactly the same mesh");
    cmp_ok($maximum_relation_error, "<=", 2e-16,
           "Phi equals (1-epsilon^2) DOS at every knot");
};

subtest "mkDOS preserves an existing table when hilb fails" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("DOS.dat", "existing table\n");
    write_file("bin/hilb", "#!/bin/sh\nexit 7\n");
    chmod(0755, "bin/hilb");

    local $ENV{PATH} = "$dir/bin:$original_path";
    isnt(system("$code/mkDOS"), 0, "mkDOS reports the backend failure");
    is(read_file("DOS.dat"), "existing table\n", "existing DOS remains intact");
    my @temporary = glob(".DOS.dat.*");
    is(scalar(@temporary), 0, "failed generation removes its temporary file");
};

subtest "Broyden remeshing matches the Steffen tool" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("source.dat", "0 0\n1 1\n2 1.5\n4 1.4\n");
    write_file(
        "mesh.dat",
        "0 0\n0.25 0\n0.75 0\n1.5 0\n2.5 0\n3.5 0\n4 0\n",
    );

    is(system("resample", "--interpolation", "steffen", "-p", "17",
              "source.dat", "mesh.dat", "expected.dat"), 0,
       "external Steffen resampling succeeds");
    my $python = <<'PYTHON';
import importlib.util
import numpy as np
import sys
spec = importlib.util.spec_from_file_location("dmft_broyden", sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
mesh = np.loadtxt(sys.argv[3], dtype=float, ndmin=2)[:, 0]
values = module.interpolate(sys.argv[2], mesh)
print("\n".join(f"{value:.17g}" for value in values))
PYTHON
    my $actual = command_output(
        "python3", "-c", $python, "$scripts/broyden.py", "source.dat", "mesh.dat"
    );
    my @actual = split(/\n/, $actual);
    my @expected = map { $_->[1] } read_table("expected.dat");
    is(scalar(@actual), scalar(@expected), "both implementations return every row");
    for my $index (0 .. $#expected) {
        cmp_ok(abs($actual[$index] - $expected[$index]), "<=", 2e-15,
               "Steffen value agrees at row " . ($index + 1));
    }
};

subtest "Bubble wrapper pins options and accepts only one scalar" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("param.loop", "clipSigma=1e-12\n");
    write_file(
        "bin/getparam",
        "#!/bin/sh\n" .
        "if [ \"\$GETPARAM_FALLBACK\" = 1 ]; then\n" .
        "  [ \"\$1\" = clip ] || exit 1\n" .
        "  printf '%s\\n' 1e-5\n" .
        "else\n" .
        "  [ \"\$1\" = clipSigma ] || exit 1\n" .
        "  printf '%s\\n' 1e-12\n" .
        "fi\n",
    );
    write_file(
        "bin/bubble",
        "#!/bin/sh\n" .
        "[ -z \"\$BUBBLE_LOG\" ] || printf '%s\\n' \"\$*\" >>\"\$BUBBLE_LOG\"\n" .
        "printf '%s\\n' \"\${BUBBLE_OUTPUT:-1.25}\"\n" .
        "exit \"\${BUBBLE_STATUS:-0}\"\n",
    );
    chmod(0755, "bin/getparam", "bin/bubble");

    local $ENV{PATH} = "$dir/bin:$original_path";
    my @options = Bubble::bubble_options(epsabs => "2e-7", epsrel => "1e-8");
    is_deeply(
        \@options,
        [
            "--gsl-error-policy", "fail",
            "--workspace-limit", "1000",
            "--interpolation", "steffen",
            "--phi-interpolation", "steffen",
            "-k", "6", "-c", "30",
            "-a", "2e-7", "-r", "1e-8", "-s", "1e-12",
        ],
        "Bubble 1.14 numerical profile is explicit",
    );
    is(Bubble::bubble_scalar("argument"), 1.25, "one finite scalar is accepted");

    {
        local $ENV{GETPARAM_FALLBACK} = 1;
        my @fallback = Bubble::bubble_options();
        cmp_ok($fallback[-1], "==", 1e-5, "legacy clip is accepted as a fallback");
    }

    local $ENV{BUBBLE_OUTPUT} = "1 2";
    my $ok = eval { Bubble::bubble_scalar("argument"); 1 };
    ok(!$ok, "multiple values are rejected");
    like($@, qr/exactly one numeric scalar/, "strict scalar error is descriptive");

    Bubble::write_scalar("value.dat", 0.125);
    is(read_file("value.dat"), "0.125\n", "scalar output is published atomically");
    my @temporary = glob("value.dat.tmp.*");
    is(scalar(@temporary), 0, "scalar publication leaves no temporary file");
};

subtest "band DOS validates before publishing Hilbert output" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin", "res");
    write_file("param.loop", "clipSigma=1e-12\n");
    write_file("mu.used", "0.5\n");
    write_file("res/resigma.dat", "0 0\n1 0\n");
    write_file("res/imsigma.dat", "0 -0.1\n1 -0.1\n");
    write_file(
        "bin/getparam",
        "#!/bin/sh\n" .
        "[ \"\$1\" = clipSigma ] || exit 1\n" .
        "printf '%s\\n' 1e-12\n",
    );
    write_file(
        "bin/hilb",
        "#!/usr/bin/env perl\n" .
        "use strict; use warnings;\n" .
        "open(my \$log, '>', \$ENV{HILB_LOG}) or die \$!;\n" .
        "print {\$log} join(' ', \@ARGV), \"\\n\"; close(\$log);\n" .
        "my (\$re, \$im) = \@ARGV[-2, -1];\n" .
        "open(my \$rf, '>', \$re) or die \$!;\n" .
        "open(my \$if, '>', \$im) or die \$!;\n" .
        "print {\$rf} \"0 -0.2\\n1 -0.3\\n\";\n" .
        "my \$spectral = \$ENV{HILB_NEGATIVE} ? -0.4 : 0.4;\n" .
        "print {\$if} \"0 \$spectral\\n1 0.5\\n\";\n" .
        "close(\$rf); close(\$if);\n",
    );
    chmod(0755, "bin/getparam", "bin/hilb");

    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{HILB_LOG} = "$dir/hilb.log";
    is(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
       "valid Hilbert output is accepted");
    like(read_file("hilb.log"), qr/-i steffen .* -x 0\.5 -c 1e-12 -f 1e-12 /,
         "band DOS pins interpolation, mu, clipping, and frequency matching");
    my $published_re = read_file("res/reaw.dat");
    my $published_im = read_file("res/imaw.dat");

    {
        local $ENV{HILB_NEGATIVE} = 1;
        isnt(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
             "negative spectral output is rejected");
    }
    is(read_file("res/reaw.dat"), $published_re, "rejected real output is not published");
    is(read_file("res/imaw.dat"), $published_im, "rejected spectral output is not published");
    my @temporary = glob("res/*.tmp.*");
    is(scalar(@temporary), 0, "rejected Hilbert output is cleaned up");
};

subtest "stable DMFT update evaluates H0 and H1 at the new mu" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin", "res");
    write_file("param.loop", "clipDelta=1e-6\nclipSigma=1e-12\n");
    write_file("mu.next", "0.25\n");
    write_file("res/resigma.dat", "0 0\n1 0\n");
    write_file("res/imsigma.dat", "0 -0.1\n1 -0.1\n");
    write_file(
        "bin/getparam",
        "#!/bin/sh\n" .
        "case \"\$1\" in\n" .
        "clipDelta) printf '%s\\n' 1e-6;;\n" .
        "clipSigma) printf '%s\\n' 1e-12;;\n" .
        "*) exit 1;;\n" .
        "esac\n",
    );
    write_file(
        "bin/hilb",
        "#!/usr/bin/env perl\n" .
        "use strict; use warnings;\n" .
        "my (\$moment, \$mu);\n" .
        "for (my \$i = 0; \$i < \@ARGV; ++\$i) {\n" .
        "  \$moment = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-n';\n" .
        "  \$mu = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-x';\n" .
        "}\n" .
        "open(my \$log, '>>', \$ENV{HILB_LOG}) or die \$!;\n" .
        "print {\$log} \"\$moment \$mu\\n\"; close(\$log);\n" .
        "my (\$re, \$im) = \@ARGV[-2, -1];\n" .
        "open(my \$rf, '>', \$re) or die \$!;\n" .
        "open(my \$if, '>', \$im) or die \$!;\n" .
        "if (\$moment == 0) { print {\$rf} \"0 1\\n1 1\\n\"; print {\$if} \"0 0\\n1 0\\n\"; }\n" .
        "else { print {\$rf} \"0 0\\n1 0\\n\"; print {\$if} \"0 -1\\n1 -1\\n\"; }\n" .
        "close(\$rf); close(\$if);\n",
    );
    chmod(0755, "bin/getparam", "bin/hilb");

    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{HILB_LOG} = "$dir/hilb.log";
    is(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
              "ReDelta.new", "ImDelta.new"), 0, "stable update succeeds");
    is(read_file("hilb.log"), "0 0.25\n1 0.25\n",
       "H0 and H1 use the same requested chemical potential");
    is(read_file("ReDelta.new"), "0 0\n1 0\n", "real ratio is F/G");
    is(read_file("ImDelta.new"), "0 -1\n1 -1\n", "imaginary ratio is causal");

    write_file("ReDelta.next.dat", "old real\n");
    write_file("ImDelta.next.dat", "old imaginary\n");
    symlink("ReDelta.next.dat", "ReDelta.dat") or die $!;
    symlink("ImDelta.next.dat", "ImDelta.dat") or die $!;
    is(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
              "ReDelta.dat", "ImDelta.dat"), 0,
       "stable update publishes through compatibility aliases");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "real Delta alias is preserved");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "imaginary Delta alias is preserved");
    is(read_file("ReDelta.next.dat"), "0 0\n1 0\n", "real alias target is updated");
    is(read_file("ImDelta.next.dat"), "0 -1\n1 -1\n", "imaginary alias target is updated");
};

subtest "zero convergence residual is valid" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("first.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("second.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("local.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("mesh.dat", "-1 0\n0 0\n1 0\n");

    is(system($^X, "$scripts/diffs", "--iteration", "2",
              "--current", "first.dat", "--previous", "second.dat",
              "--local", "local.dat", "--mesh", "mesh.dat",
              "--output-dir", "."), 0, "diffs accepts an exact zero residual");
    like(read_file("DIFFS_C"), qr/^2\s+0(?:\.0*)?\s*$/,
         "zero consecutive residual is recorded");
    like(read_file("DIFFS_LatLoc"), qr/^2\s+0(?:\.0*)?\s*$/,
         "zero lattice-local residual is recorded");
    ok(!-e ".diff-residual.tmp", "temporary residual is removed");
};

subtest "KK and optical sum-rule utilities use strict Steffen tools" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("spectra");
    my $spectral = "-1 0\n-0.5 1\n0.5 1\n1 0\n";
    write_file("spectra/c-imF.dat", $spectral);
    write_file("spectra/c-imG.dat", $spectral);
    write_file("spectra/c-imI.dat", $spectral);
    is(system($^X, "$scripts/realparts", "spectra"), 0,
       "all three real parts use the installed Steffen KK backend");
    is(scalar(read_table("spectra/c-reG.dat")), 4, "KK preserves the input mesh");

    write_file("ekin.dat", "-1\n");
    write_file("cond.opt-PHI.dat", "0 1\n1 1\n2 1\n");
    my $sum_rule = command_output($^X, "$scripts/sumrule");
    $sum_rule =~ /\Ar=([^\s]+)\s*\z/ or die "Unexpected sumrule output: $sum_rule";
    my $expected = 8 / (3 * (4 * atan2(1, 1)) ** 2);
    cmp_ok(abs($1 - $expected), "<=", 2e-15,
           "sumrule integrates the represented optical curve");
};

subtest "asymmetric first-moment publication stops on Bubble failure" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("param.eps-bubble", "existing\n");
    write_file("bin/bubble", "#!/bin/sh\nprintf partial\nexit 9\n");
    chmod(0755, "bin/bubble");

    local $ENV{PATH} = "$dir/bin:$original_path";
    isnt(system("bash", "$code/../more_examples/asym/mkparameps-bubble"), 0,
         "Bubble failure propagates from the wrapper");
    is(read_file("param.eps-bubble"), "existing\n", "existing moment remains intact");
    my @temporary = glob("param.eps-bubble.tmp.*");
    is(scalar(@temporary), 0, "failed moment output is removed");
};

subtest "directory comparison uses the Steffen convergence norm" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("reference");
    my $table = "-1 0.1\n0 0.2\n1 0.1\n";
    write_file("spectrum.dat", $table);
    write_file("reference/spectrum.dat", $table);
    write_file("scalar.dat", "1.25\n");
    write_file("reference/scalar.dat", "1.25\n");

    is(system($^X, "$code/../extra/compare_dir", "reference"), 0,
       "identical scalar and table files compare successfully");
    write_file("scalar.dat", "2.25\n");
    isnt(system($^X, "$code/../extra/compare_dir", "reference"), 0,
         "a negative signed scalar difference cannot cancel the error norm");
};

subtest "installed Hilbert and Bubble backends complete strict profiles" => sub {
    plan skip_all => "hilb is not installed" unless executable_available("hilb");
    plan skip_all => "bubble is not installed" unless executable_available("bubble");

    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    symlink("$code/DOS.dat", "DOS.dat") or die $!;
    symlink("$code/PHI.dat", "PHI.dat") or die $!;
    write_file("param.loop", "clipDelta=1e-6\nclipSigma=1e-12\n");
    write_file("param.mu", "0\n");
    my ($real_sigma, $imaginary_sigma) = ("", "");
    for my $index (-200 .. 200) {
        my $omega = $index / 100;
        $real_sigma .= sprintf("%.17g 0\n", $omega);
        $imaginary_sigma .= sprintf("%.17g -0.05\n", $omega);
    }
    write_file("resigma.dat", $real_sigma);
    write_file("imsigma.dat", $imaginary_sigma);

    is(system($^X, "$scripts/bandDOS", ".", "param.mu"), 0,
       "real analytic Hilbert backend produces a causal lattice spectrum");
    is(system($^X, "$scripts/dmftDOS-stable", ".", "param.mu",
              "ReDelta.dat", "ImDelta.dat"), 0,
       "real H0/H1 transforms produce a stable hybridization");
    my @spectrum = read_table("imaw.dat");
    my @hybridization = read_table("ImDelta.dat");
    is(scalar(@spectrum), 401, "Hilbert output preserves the self-energy mesh");
    ok(!grep({ $_->[1] < 0 } @spectrum), "lattice spectrum is nonnegative");
    ok(!grep({ $_->[1] >= 0 } @hybridization), "hybridization remains causal");

    my @standard_options = Bubble::bubble_options(
        epsabs => "1e-9", epsrel => "1e-8"
    );
    my $standard_dc = Bubble::bubble_scalar(
        @standard_options, 5, 2, 0, 0.02, 0,
        "resigma.dat", "imsigma.dat",
    );
    cmp_ok($standard_dc, ">", 0, "standard Bethe DC profile completes");

    my @dc_options = Bubble::bubble_options(epsabs => "1e-7", epsrel => "1e-8");
    my $dc = Bubble::bubble_scalar(
        @dc_options, "-p", "PHI.dat", 0, 2, 0, 0.02, 0,
        "resigma.dat", "imsigma.dat",
    );
    cmp_ok($dc, ">", 0, "tabulated-Phi DC profile completes");

    my @kinetic_options = Bubble::bubble_options(
        epsabs => "1e-7", epsrel => "1e-8"
    );
    my $kinetic = Bubble::bubble_scalar(
        @kinetic_options, "-e", "1", "-p", "DOS.dat", "-f",
        0, 1, 0, 0.02, 0, "resigma.dat", "imsigma.dat",
    );
    cmp_ok($kinetic, "<", 0, "occupied kinetic-energy profile completes");

    my @optical_options = Bubble::bubble_options(
        epsabs => "2e-7", epsrel => "1e-8"
    );
    my $optical = Bubble::bubble_scalar(
        "-O", "0.1", @optical_options, "-p", "PHI.dat",
        0, 2, 0, 0.02, 0, "resigma.dat", "imsigma.dat",
    );
    cmp_ok($optical, ">", 0, "positive-frequency optical profile completes");
};

chdir($original) or die $!;
done_testing();

sub max {
    return $_[0] > $_[1] ? $_[0] : $_[1];
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
        @fields == 2 or die "Expected two columns in $path at line $.";
        push(@rows, [map { 0.0 + $_ } @fields]);
    }
    close($fh) or die "Can't close $path: $!";
    return @rows;
}

sub command_output {
    my @command = @_;
    open(my $fh, "-|", @command) or die "Can't run $command[0]: $!";
    local $/;
    my $output = <$fh>;
    close($fh) or die "$command[0] failed: $?";
    return defined($output) ? $output : "";
}

sub executable_available {
    my ($name) = @_;
    for my $directory (split(/:/, $ENV{PATH})) {
        return 1 if -x "$directory/$name";
    }
    return 0;
}
