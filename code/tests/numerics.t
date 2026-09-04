#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use Digest::SHA qw(sha256_hex);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use Test::More;

my $code = "$FindBin::Bin/..";
my $scripts = "$code/scripts";
my $original = getcwd();
my $original_path = $ENV{PATH};
my $external_mode = $ENV{DMFT_TEST_EXTERNAL} // "auto";
$external_mode =~ /\A(?:0|1|auto)\z/
    or die "DMFT_TEST_EXTERNAL must be 0, 1, or auto\n";

require "$scripts/Bubble.pm";

subtest "Bethe tables use one normalized edge-refined mesh" => sub {
    require_external_tools("hilb") or return;
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
    like(read_file("$code/param.loop"), qr/^density_interpolation=steffen$/m,
         "solver configuration pins Steffen density interpolation");
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

subtest "mkDOS reports warnings and preserves output on failure" => sub {
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

    write_file(
        "bin/hilb",
        "#!/bin/sh\n" .
        "printf '%s\\n' 'algorithm=analytic-piecewise-polynomial' " .
        "'dos.integral=1' 'WARNING - qag error: 18 -- roundoff error' >&2\n",
    );
    chmod(0755, "bin/hilb");
    open(my $saved_stderr, ">&", \*STDERR) or die $!;
    open(STDERR, ">", "mkdos.stderr") or die $!;
    my $status = system("$code/mkDOS");
    open(STDERR, ">&", $saved_stderr) or die $!;
    close($saved_stderr);
    is($status, 0, "mkDOS accepts a finite warning-mode normalization");
    like(read_file("mkdos.stderr"), qr/WARNING - qag error: 18/,
         "successful hilb warning remains visible");
    is(scalar(read_table("DOS.dat")), 2601, "warning-mode DOS is published");
};

subtest "Broyden remeshing matches the Steffen tool" => sub {
    require_external_tools("resample") or return;
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

subtest "causal Delta is projected and reconstructed from Gamma" => sub {
    require_external_tools(qw(getparam kk)) or return;
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("param.loop", "clipDelta=0.1\n");
    write_file("param.eps", "0.25\n");
    write_file(
        "Gamma.raw.dat",
        "-2 1e-9\n-1 0.05\n-0.5 0.2\n0.5 -5e-8\n1 0.4\n2 2e-9\n",
    );

    is(system($^X, "$scripts/causalDelta", "Gamma.raw.dat", "Delta.dat",
              "ReDelta.dat", "ImDelta.dat"), 0,
       "causal reconstruction succeeds");
    my @gamma = read_table("Delta.dat");
    my @imaginary = read_table("ImDelta.dat");
    my @real = read_table("ReDelta.dat");
    is_deeply([map { $_->[1] } @gamma], [0, 0.1, 0.2, 0.1, 0.4, 0],
              "only tolerance-bounded Gamma values are projected");
    for my $index (0 .. $#gamma) {
        is($imaginary[$index][0], $gamma[$index][0],
           "imaginary mesh agrees at row " . ($index + 1));
        cmp_ok(abs($imaginary[$index][1] + $gamma[$index][1]), "<=", 1e-16,
               "ImDelta is exactly -Gamma at row " . ($index + 1));
    }

    is(system("kk", "--interpolation", "steffen", "ImDelta.dat", "dynamic.dat"),
       0, "independent Steffen KK succeeds");
    my @dynamic = read_table("dynamic.dat");
    for my $index (0 .. $#real) {
        cmp_ok(abs($real[$index][1] - $dynamic[$index][1] - 0.25), "<=", 2e-15,
               "ReDelta is KK[ImDelta]+param.eps at row " . ($index + 1));
    }

    my $published_gamma = read_file("Delta.dat");
    my $published_real = read_file("ReDelta.dat");
    my $published_imaginary = read_file("ImDelta.dat");
    write_file("Gamma.raw.dat", "-1 0\n0 broken\n1 0\n");
    isnt(system($^X, "$scripts/causalDelta", "Gamma.raw.dat", "Delta.dat",
                "ReDelta.dat", "ImDelta.dat"), 0,
         "malformed Gamma is rejected");
    is(read_file("Delta.dat"), $published_gamma,
       "failed reconstruction preserves Gamma");
    is(read_file("ReDelta.dat"), $published_real,
       "failed reconstruction preserves ReDelta");
    is(read_file("ImDelta.dat"), $published_imaginary,
       "failed reconstruction preserves ImDelta");
    write_file(
        "Gamma.raw.dat",
        "-2 0\n-1 0.2\n-0.5 0.2\n0.5 0.2\n1 0.2\n2 0\n",
    );
    write_file("param.eps", "1e999\n");
    isnt(system($^X, "$scripts/causalDelta", "Gamma.raw.dat", "Delta.dat",
                "ReDelta.dat", "ImDelta.dat"), 0,
         "non-finite param.eps is rejected");
    is(read_file("Delta.dat"), $published_gamma,
       "invalid param.eps preserves the published triplet");
    my @temporary = glob("*.tmp.*");
    is(scalar(@temporary), 0, "causal reconstruction leaves no temporary files");
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
            "--gsl-error-policy", "warn",
            "--workspace-limit", "1000",
            "--interpolation", "steffen",
            "--phi-interpolation", "steffen",
            "-k", "6", "-c", "30",
            "-a", "2e-7", "-r", "1e-8", "-s", "1e-12",
        ],
        "Bubble 1.14 warning-mode numerical profile is explicit",
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

    {
        local $ENV{BUBBLE_OUTPUT} = "1e999";
        my $finite = eval { Bubble::bubble_scalar("argument"); 1 };
        ok(!$finite, "overflow-form Bubble output is rejected as non-finite");
    }

    Bubble::write_scalar("value.dat", 0.125);
    is(read_file("value.dat"), "0.125\n", "scalar output is published atomically");
    my @temporary = glob("value.dat.tmp.*");
    is(scalar(@temporary), 0, "scalar publication leaves no temporary file");
};

subtest "lattice DOS validates warning-mode Bubble tables" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("param.loop", "T=0.02\nclipSigma=1e-12\n");
    write_file("param.mu", "0\n");
    write_file("resigma.dat", "0 0\n1 0\n");
    write_file("imsigma.dat", "0 -0.1\n1 -0.1\n");
    write_file(
        "bin/getparam",
        "#!/bin/sh\n" .
        "case \"\$1\" in\n" .
        "T) printf '%s\\n' 0.02;;\n" .
        "clipSigma) printf '%s\\n' 1e-12;;\n" .
        "*) exit 1;;\n" .
        "esac\n",
    );
    write_file(
        "bin/bubble",
        "#!/usr/bin/env perl\n" .
        "use strict; use warnings;\n" .
        "join(' ', \@ARGV) =~ /--gsl-error-policy warn/ or exit 91;\n" .
        "my \$output;\n" .
        "for (my \$i = 0; \$i < \@ARGV; ++\$i) {\n" .
        "  \$output = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-o';\n" .
        "}\n" .
        "open(my \$fh, '>', \$output) or die \$!;\n" .
        "if (defined \$ENV{BUBBLE_TABLE}) { print {\$fh} \$ENV{BUBBLE_TABLE}; }\n" .
        "else {\n" .
        "  for my \$i (0 .. 10000) {\n" .
        "    my \$omega = \$ENV{BUBBLE_WRONG_MESH} ? \$i : \$i / 10000;\n" .
        "    printf {\$fh} \"%.17g 0.5\\n\", \$omega;\n" .
        "  }\n" .
        "}\n" .
        "close(\$fh);\n" .
        "print STDERR \"WARNING - qag error: 18 -- roundoff error\\n\";\n",
    );
    chmod(0755, "bin/getparam", "bin/bubble");

    local $ENV{PATH} = "$dir/bin:$original_path";
    is(system($^X, "$scripts/dos"), 0,
       "finite Bubble table is published despite a visible warning");
    my $published = read_file("ldos.dat");

    {
        local $ENV{BUBBLE_WRONG_MESH} = 1;
        isnt(system($^X, "$scripts/dos"), 0,
             "Bubble table on the wrong interval remains fatal");
    }
    is(read_file("ldos.dat"), $published,
       "wrong-interval Bubble table is not published");

    {
        local $ENV{BUBBLE_TABLE} = "0 0.5\n1 0.4\n";
        isnt(system($^X, "$scripts/dos"), 0,
             "truncated Bubble table remains fatal");
    }
    is(read_file("ldos.dat"), $published, "truncated Bubble table is not published");

    {
        local $ENV{BUBBLE_TABLE} = "0 1e999\n1 0.4\n";
        isnt(system($^X, "$scripts/dos"), 0,
             "non-finite Bubble table remains fatal");
    }
    is(read_file("ldos.dat"), $published, "invalid Bubble table is not published");
};

subtest "band DOS validates before publishing Hilbert output" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin", "res");
    write_file("param.loop", "clipSigma=1e-12\n");
    write_file("mu.used", "0.5\n");
    write_file("DOS.dat", "0 1\n1 1\n");
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
        "my \$real = \$ENV{HILB_TRUNCATED} ? \"0 -0.2\\n\" : \"0 -0.2\\n1 -0.3\\n\";\n" .
        "print {\$rf} \$real;\n" .
        "my \$spectral = \$ENV{HILB_NEGATIVE} ? -0.4 : 0.4;\n" .
        "my \$imaginary = \$ENV{HILB_TRUNCATED} " .
        "? \"0 \$spectral\\n\" : \"0 \$spectral\\n1 0.5\\n\";\n" .
        "print {\$if} \$imaginary;\n" .
        "close(\$rf); close(\$if);\n",
    );
    chmod(0755, "bin/getparam", "bin/hilb");

    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{HILB_LOG} = "$dir/hilb.log";
    is(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
       "valid Hilbert output is accepted");
    like(read_file("hilb.log"),
         qr/-i steffen --gsl-error-policy warn .* -x 0\.5 -c 1e-12 -f 1e-12 /,
         "band DOS pins warning policy, interpolation, mu, and clipping");
    my $published_re = read_file("res/reaw.dat");
    my $published_im = read_file("res/imaw.dat");

    is(system(
        $^X, "$scripts/bandDOS", "--mu-value", "0.25",
        "--re-output", "trial.re", "--im-output", "trial.im",
        "--dos", "DOS.dat", "res",
    ), 0, "direct mu and explicit output paths are accepted");
    like(read_file("hilb.log"), qr/-d DOS\.dat -x 0\.25 /,
         "trial evaluation forwards the selected DOS and mu");
    is(read_file("trial.re"), $published_re, "explicit real output is published");
    is(read_file("trial.im"), $published_im, "explicit spectral output is published");

    write_file("real.target", "old real\n");
    write_file("imaginary.target", "old imaginary\n");
    symlink("real.target", "real.alias") or die $!;
    symlink("imaginary.target", "imaginary.alias") or die $!;
    is(system(
        $^X, "$scripts/bandDOS", "--mu-value", "0.25",
        "--re-output", "real.alias", "--im-output", "imaginary.alias", "res",
    ), 0, "explicit outputs publish through compatibility symlinks");
    is(readlink("real.alias"), "real.target", "real output symlink is preserved");
    is(readlink("imaginary.alias"), "imaginary.target",
       "spectral output symlink is preserved");

    my $sigma_before = read_file("res/resigma.dat");
    isnt(system(
        $^X, "$scripts/bandDOS", "--mu-value", "0.25",
        "--re-output", "res/resigma.dat", "--im-output", "collision.im", "res",
    ), 0, "trial output cannot overwrite a frozen-Sigma input");
    is(read_file("res/resigma.dat"), $sigma_before,
       "output collision preserves the self-energy input");

    {
        local $ENV{HILB_NEGATIVE} = 1;
        isnt(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
             "negative spectral output is rejected");
    }
    is(read_file("res/reaw.dat"), $published_re, "rejected real output is not published");
    is(read_file("res/imaw.dat"), $published_im, "rejected spectral output is not published");

    {
        local $ENV{HILB_TRUNCATED} = 1;
        isnt(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
             "truncated warning-mode Hilbert output is rejected");
    }
    is(read_file("res/reaw.dat"), $published_re, "truncated real output is not published");
    is(read_file("res/imaw.dat"), $published_im,
       "truncated spectral output is not published");
    write_file("mu.used", "0.5\njunk\n");
    isnt(system($^X, "$scripts/bandDOS", "res", "mu.used"), 0,
         "corrupt chemical-potential input is rejected");
    is(read_file("res/imaw.dat"), $published_im,
       "corrupt chemical potential does not replace output");
    my @temporary = glob("res/*.tmp.*");
    is(scalar(@temporary), 0, "rejected Hilbert output is cleaned up");
};

subtest "stable DMFT update evaluates H0 and H1 at the new mu" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin", "res");
    write_file("param.loop", "clipDelta=1e-6\nclipSigma=1e-12\n");
    write_file("mu.next", "0.25\n");
    write_file("DOS.dat", "0 1\n1 1\n");
    write_file("res/resigma.dat", "0 0\n1 0\n");
    write_file("res/imsigma.dat", "0 -0.1\n1 -0.1\n");
    write_file(
        "bin/getparam",
        "#!/bin/sh\n" .
        "case \"\$1\" in\n" .
        "clipDelta) printf '%s\\n' \"\${CLIP_DELTA:-1e-6}\";;\n" .
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
        "if (\$ENV{HILB_EMPTY}) { close(\$rf); close(\$if); exit 0; }\n" .
        "if (\$moment == 0) { print {\$rf} \"0 1\\n1 1\\n\"; print {\$if} \"0 0\\n1 0\\n\"; }\n" .
        "else { print {\$rf} \"0 0\\n1 0\\n\"; print {\$if} \"0 -1\\n1 -1\\n\"; }\n" .
        "close(\$rf); close(\$if);\n",
    );
    chmod(0755, "bin/getparam", "bin/hilb");

    local $ENV{PATH} = "$dir/bin:$original_path";
    local $ENV{HILB_LOG} = "$dir/hilb.log";
    is(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
              "Delta.new"), 0, "stable update succeeds");
    is(read_file("hilb.log"), "0 0.25\n1 0.25\n",
       "H0 and H1 use the same requested chemical potential");
    is(read_file("Delta.new"), "0 1\n1 1\n", "raw update stores Gamma=-Im(F/G)");

    write_file("H0.re.dat", "0 1\n1 1\n");
    write_file("H0.im.dat", "0 0\n1 0\n");
    write_h0_metadata("H0.mu", 0.25, 1e-12, "res/resigma.dat",
                      "res/imsigma.dat", "DOS.dat", "H0.re.dat", "H0.im.dat");
    unlink("hilb.log") or die $!;
    is(system(
        $^X, "$scripts/dmftDOS-stable",
        "--h0-real", "H0.re.dat", "--h0-imaginary", "H0.im.dat",
        "--h0-mu", "H0.mu", "--dos", "DOS.dat",
        "res", "mu.next", "Delta.cached",
    ), 0, "stable update accepts an accurate occupancy H0 cache");
    is(read_file("hilb.log"), "1 0.25\n",
       "cached update skips H0 but still evaluates H1 at the same mu");
    is(read_file("Delta.cached"), "0 1\n1 1\n",
       "cached and uncached raw updates agree");

    write_h0_metadata("H0.mu", 0.5, 1e-12, "res/resigma.dat",
                      "res/imsigma.dat", "DOS.dat", "H0.re.dat", "H0.im.dat");
    write_file("Delta.cached", "published cache result\n");
    isnt(system(
        $^X, "$scripts/dmftDOS-stable",
        "--h0-real", "H0.re.dat", "--h0-imaginary", "H0.im.dat",
        "--h0-mu", "H0.mu", "res", "mu.next", "Delta.cached",
    ), 0, "stale H0 cache is rejected");
    is(read_file("Delta.cached"), "published cache result\n",
       "stale cache does not replace raw Gamma");
    write_h0_metadata("H0.mu", 0.25, 1e-12, "res/resigma.dat",
                      "res/imsigma.dat", "DOS.dat", "H0.re.dat", "H0.im.dat");
    write_file("res/resigma.dat", "0 0.25\n1 0.25\n");
    isnt(system(
        $^X, "$scripts/dmftDOS-stable",
        "--h0-real", "H0.re.dat", "--h0-imaginary", "H0.im.dat",
        "--h0-mu", "H0.mu", "res", "mu.next", "Delta.cached",
    ), 0, "H0 cache from a different frozen self-energy is rejected");
    is(read_file("Delta.cached"), "published cache result\n",
       "cache provenance failure preserves raw Gamma");
    write_file("res/resigma.dat", "0 0\n1 0\n");

    write_file("Delta.next.dat", "old Gamma\n");
    symlink("Delta.next.dat", "Delta.dat") or die $!;
    is(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
              "Delta.dat"), 0,
       "stable update publishes through compatibility aliases");
    is(readlink("Delta.dat"), "Delta.next.dat", "Gamma alias is preserved");
    is(read_file("Delta.next.dat"), "0 1\n1 1\n", "Gamma alias target is updated");

    {
        local $ENV{HILB_EMPTY} = 1;
        isnt(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
                    "Delta.dat"), 0,
             "empty warning-mode Hilbert output is rejected");
    }
    is(read_file("Delta.next.dat"), "0 1\n1 1\n",
       "empty output does not replace Gamma");
    write_file("mu.next", "0.25\njunk\n");
    isnt(system($^X, "$scripts/dmftDOS-stable", "res", "mu.next",
                "Delta.dat"), 0,
          "corrupt chemical-potential input is rejected");
    is(read_file("Delta.next.dat"), "0 1\n1 1\n",
       "corrupt chemical potential does not replace Gamma");
};

subtest "warning-mode convergence estimates continue" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("first.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("second.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("local.dat", "-1 0.1\n0 0.2\n1 0.1\n");
    write_file("mesh.dat", "-1 0\n0 0\n1 0\n");
    write_file(
        "bin/integ",
        "#!/bin/sh\n" .
        "case \" \$* \" in\n" .
        "  *\" --gsl-error-policy warn \"*) ;;\n" .
        "  *) printf '%s\\n' 'missing warning policy' >&2; exit 91;;\n" .
        "esac\n" .
        "printf '%s\\n' 'WARNING - qag error: 18 -- roundoff error' >&2\n" .
        "printf '%s\\n' \"\${INTEG_OUTPUT:-0}\"\n",
    );
    chmod(0755, "bin/integ");
    local $ENV{PATH} = "$dir/bin:$original_path";

    is(system($^X, "$scripts/diffs", "--iteration", "2",
              "--current", "first.dat", "--previous", "second.dat",
              "--local", "local.dat", "--mesh", "mesh.dat",
              "--output-dir", "."), 0,
       "diffs records finite best estimates accompanied by warnings");
    like(read_file("DIFFS_C"), qr/^2\s+0(?:\.0*)?\s*$/,
         "zero consecutive residual is recorded");
    like(read_file("DIFFS_LatLoc"), qr/^2\s+0(?:\.0*)?\s*$/,
         "zero lattice-local residual is recorded");
    ok(!-e ".diff-residual.tmp", "temporary residual is removed");

    {
        local $ENV{INTEG_OUTPUT} = "1e999";
        isnt(system($^X, "$scripts/diffs", "--iteration", "2",
                    "--current", "first.dat", "--previous", "second.dat",
                    "--local", "local.dat", "--mesh", "mesh.dat",
                    "--output-dir", "."), 0,
             "non-finite best estimates remain fatal");
    }
};

subtest "Delta support check removes only the represented floor" => sub {
    require_external_tools(qw(getparam integ)) or return;
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("param.loop", "bandrescale=0.5\nclipDelta=0.1\n");
    write_file("Delta.dat", "-2 0\n-1 0.1\n1 0.1\n2 0\n");
    is(system($^X, "$scripts/checkDelta"), 0,
       "pure numerical floor outside the solver window is ignored");
    my @floor_temporary = glob(".Delta.floor.tmp.*");
    is(scalar(@floor_temporary), 0, "floor integration table is removed");

    write_file("Delta.dat", "-2 0\n-1 0.2\n1 0.2\n2 0\n");
    isnt(system($^X, "$scripts/checkDelta"), 0,
         "physical hybridization weight outside the window remains fatal");
};

subtest "KK and optical sum-rule utilities use configured Steffen tools" => sub {
    require_external_tools(qw(kk integ)) or return;
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

    write_file("ekin.dat", "1e999\n");
    my $finite = eval { command_output($^X, "$scripts/sumrule"); 1 };
    ok(!$finite, "sumrule rejects overflow-form scalar input");
};

subtest "asymmetric first-moment publication stops on subprocess failure" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("param.eps-bubble", "existing\n");
    write_file(
        "bin/bubble",
        "#!/bin/sh\n" .
        "case \" \$* \" in\n" .
        "  *\" --gsl-error-policy warn \"*) ;;\n" .
        "  *) exit 91;;\n" .
        "esac\n" .
        "[ -z \"\$BUBBLE_WARNING\" ] || printf '%s\\n' " .
        "'WARNING - qag error: 18 -- roundoff error' >&2\n" .
        "printf '%s\\n' \"\${BUBBLE_OUTPUT:-partial}\"\n" .
        "exit \"\${BUBBLE_STATUS:-9}\"\n",
    );
    chmod(0755, "bin/bubble");

    local $ENV{PATH} = "$dir/bin:$original_path";
    isnt(system("bash", "$code/../more_examples/asym/mkparameps-bubble"), 0,
         "Bubble failure propagates from the wrapper");
    is(read_file("param.eps-bubble"), "existing\n", "existing moment remains intact");
    my @temporary = glob("param.eps-bubble.tmp.*");
    is(scalar(@temporary), 0, "failed moment output is removed");

    {
        local $ENV{BUBBLE_STATUS} = 0;
        local $ENV{BUBBLE_OUTPUT} = "not-a-number";
        isnt(system("bash", "$code/../more_examples/asym/mkparameps-bubble"), 0,
             "malformed warning-mode output remains fatal");
    }
    is(read_file("param.eps-bubble"), "existing\n",
       "malformed output is not published");

    {
        local $ENV{BUBBLE_STATUS} = 0;
        local $ENV{BUBBLE_OUTPUT} = "0.125";
        local $ENV{BUBBLE_WARNING} = 1;
        is(system("bash", "$code/../more_examples/asym/mkparameps-bubble"), 0,
           "finite best estimate is published despite a visible warning");
    }
    is(read_file("param.eps-bubble"), "0.125\n",
       "finite warning-mode result is published");
};

subtest "directory comparison uses the Steffen convergence norm" => sub {
    require_external_tools(qw(resample subtracty integ subtract)) or return;
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

subtest "installed Hilbert and Bubble backends complete warning-mode profiles" => sub {
    require_external_tools(qw(getparam hilb kk bubble)) or return;

    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    symlink("$code/DOS.dat", "DOS.dat") or die $!;
    symlink("$code/PHI.dat", "PHI.dat") or die $!;
    write_file("param.loop", "clipDelta=1e-6\nclipSigma=1e-12\n");
    write_file("param.eps", "0\n");
    write_file("param.mu", "0\n");
    my ($real_sigma, $imaginary_sigma) = ("", "");
    for my $index (-200 .. -1, 1 .. 200) {
        my $omega = $index / 100;
        $real_sigma .= sprintf("%.17g 0\n", $omega);
        $imaginary_sigma .= sprintf("%.17g -0.05\n", $omega);
    }
    write_file("resigma.dat", $real_sigma);
    write_file("imsigma.dat", $imaginary_sigma);

    is(system($^X, "$scripts/bandDOS", ".", "param.mu"), 0,
       "real analytic Hilbert backend produces a causal lattice spectrum");
    is(system($^X, "$scripts/dmftDOS-stable", ".", "param.mu",
              "Delta.raw.dat"), 0,
       "real H0/H1 transforms produce a raw Gamma update");
    is(system($^X, "$scripts/causalDelta", "Delta.raw.dat", "Delta.dat",
              "ReDelta.dat", "ImDelta.dat"), 0,
       "raw Gamma is projected into a causal hybridization");
    my @spectrum = read_table("imaw.dat");
    my @gamma = read_table("Delta.dat");
    my @hybridization = read_table("ImDelta.dat");
    is(scalar(@spectrum), 400, "Hilbert output preserves the self-energy mesh");
    is(scalar(@hybridization), 400,
       "stable hybridization preserves the self-energy mesh");
    ok(!grep({ $_->[1] < 0 } @spectrum), "lattice spectrum is nonnegative");
    is($gamma[0][1], 0, "lower Gamma support guard is zero");
    is($gamma[-1][1], 0, "upper Gamma support guard is zero");
    ok(!grep({ $_->[1] < 1e-6 } @gamma[1 .. $#gamma - 1]),
       "interior Gamma respects clipDelta");
    my $consistent_imaginary = 1;
    for my $index (0 .. $#hybridization) {
        $consistent_imaginary = 0
            if abs($hybridization[$index][1] + $gamma[$index][1]) > 1e-14;
    }
    ok($consistent_imaginary, "imaginary hybridization equals -Gamma");

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

sub require_external_tools {
    my @tools = @_;
    if ($external_mode eq "0") {
        plan skip_all => "external numerical tools disabled by DMFT_TEST_EXTERNAL=0";
        return 0;
    }

    my @missing = grep { !executable_available($_) } @tools;
    return 1 if !@missing;

    my $message = "missing external numerical tools: " . join(", ", @missing);
    die "$message\n" if $external_mode eq "1";
    plan skip_all => $message;
    return 0;
}

sub write_h0_metadata {
    my ($path, $mu, $clip, $resigma, $imsigma, $dos, $h0_real, $h0_imaginary) = @_;
    write_file(
        $path,
        "version=1\n" .
        "mu=$mu\n" .
        "clip_sigma=$clip\n" .
        "resigma_sha256=" . sha256_hex(read_file($resigma)) . "\n" .
        "imsigma_sha256=" . sha256_hex(read_file($imsigma)) . "\n" .
        "dos_sha256=" . sha256_hex(read_file($dos)) . "\n" .
        "h0_real_sha256=" . sha256_hex(read_file($h0_real)) . "\n" .
        "h0_imaginary_sha256=" . sha256_hex(read_file($h0_imaginary)) . "\n",
    );
}
