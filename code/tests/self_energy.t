#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use List::Util qw(max);
use Test::More;

my $sigmatrick = "$FindBin::Bin/../scripts/sigmatrick";
my $original = getcwd();
my $original_path = $ENV{PATH};
my $pi = 4 * atan2(1, 1);
my @omega = (-1, -0.25, 0.25, 1);
my @gamma = (0, 0.1, 0.1, 0);
my @outputs = qw(c-self.dat imsigma.dat resigma.dat);
my $external_mode = $ENV{DMFT_TEST_EXTERNAL} // "auto";
$external_mode =~ /\A(?:0|1|auto)\z/
    or die "DMFT_TEST_EXTERNAL must be 0, 1, or auto\n";

subtest "improved estimator and Dyson reconstruction" => sub {
    my $dir = self_energy_fixture(clip => 1e-6, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(
        G => [1, -1],
        F => [0.5, -0.25],
        I => [0.2, -0.3],
    );
    stage_production_layout();
    write_file("res/param.loop", "sigmatrick=old\n");

    is(run_sigmatrick("res"), 0, "synthetic estimator succeeds in the staged layout");
    my @real = read_table("res/resigma.dat");
    my @imaginary = read_table("res/imsigma.dat");
    my @spectral = read_table("res/c-self.dat");
    is(readlink("res/ReDelta.dat"), "ReDelta.used.dat",
       "production ReDelta alias is preserved");
    is(scalar(@real), scalar(@omega), "self-energy covers the full mesh");

    my $expected_real = 0.08125;
    my $expected_imaginary = -0.26875;
    for my $index (0 .. $#omega) {
        is($real[$index][0], $omega[$index], "real mesh row " . ($index + 1));
        is($imaginary[$index][0], $omega[$index],
           "imaginary mesh row " . ($index + 1));
        cmp_ok(abs($real[$index][1] - $expected_real), "<=", 2e-14,
               "estimator real part row " . ($index + 1));
        cmp_ok(abs($imaginary[$index][1] - $expected_imaginary), "<=", 2e-14,
               "estimator imaginary part row " . ($index + 1));

        my $denominator_real = $omega[$index] + 0.2 - 0.05 - $expected_real;
        my $denominator_imaginary = $gamma[$index] - $expected_imaginary;
        my $expected_spectral = $denominator_imaginary / (
            $pi * ($denominator_real**2 + $denominator_imaginary**2)
        );
        cmp_ok(abs($spectral[$index][1] - $expected_spectral), "<=", 2e-14,
               "Dyson spectrum row " . ($index + 1));
    }

    my @default = map { read_file("res/$_") } @outputs;
    write_file("param.loop", "clipSigma=1e-6\nsigmatrick=new\n" .
               "w_crossover=invalid\nw_width=invalid\n");
    is(run_sigmatrick("res"), 0, "explicit new reads working-directory parameters");
    is_deeply([map { read_file("res/$_") } @outputs], \@default,
              "absent mode and explicit new are byte-identical");
    write_file("param.loop", "clip=1e-6\n  sigmatrick = NEW  \n");
    is(run_sigmatrick("res"), 0, "legacy clip fallback and whitespace are accepted");
    is_deeply([map { read_file("res/$_") } @outputs], \@default,
              "clip fallback and case-normalized mode preserve the result");
    ok(!-e "kk.args", "new does not invoke KK or validate crossover parameters");
};

subtest "positive imaginary self-energy is clipped" => sub {
    my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(
        G => [1, -1],
        F => [0, 0],
        I => [0, 0.2],
    );

    is(run_sigmatrick(), 0, "causality clipping succeeds");
    my @real = read_table("resigma.dat");
    my @imaginary = read_table("imsigma.dat");
    my @spectral = read_table("c-self.dat");
    for my $index (0 .. $#omega) {
        cmp_ok(abs($real[$index][1] - 0.1), "<=", 2e-14,
               "clipping preserves the real part at row " . ($index + 1));
        cmp_ok(abs($imaginary[$index][1] + 0.01), "<=", 2e-14,
               "imaginary part is clipped at row " . ($index + 1));
        cmp_ok($spectral[$index][1], ">", 0,
               "clipped Dyson spectrum is positive at row " . ($index + 1));
    }
};

subtest "old is F/G without I, Hartree, or KK" => sub {
    my $dir = self_energy_fixture(clip => 0.01, sigma_h => 4, mode => "old",
                                 parameters => "w_crossover=bad\nw_width=bad\n");
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(G => [1, -1], F => [0.25, -0.75], I => [8, -9]);
    stage_production_layout();
    write_file("res/param.loop", "sigmatrick=new\n");

    is(run_sigmatrick("res"), 0, "old succeeds with a nonzero Hartree input");
    assert_result([(0.5) x @omega], [(-0.25) x @omega], "old", "res/");
    my @published = map { read_file("res/$_") } @outputs;
    unlink("res/c-reI.dat", "res/c-imI.dat", "res/customfdm.avg",
           "bin/extractcolumn") == 4 or die $!;
    is(run_sigmatrick("res"), 0, "old needs neither I nor Hartree extraction");
    is_deeply([map { read_file("res/$_") } @outputs], \@published,
              "old is independent of I and does not double-count Hartree");
    write_table("res/c-reF.dat", \@omega, [(-0.75 / $pi) x @omega]);
    write_table("res/c-imF.dat", \@omega, [(0.25 / $pi) x @omega]);
    is(run_sigmatrick("res"), 0, "old clips a positive imaginary F/G");
    assert_result([(0.5) x @omega], [(-0.01) x @omega], "clipped old", "res/");
    ok(!-e "kk.args", "old does not invoke KK");
};

subtest "new-im uses normalized imaginary correlators and clips before KK" => sub {
    my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.4, mode => "new-im",
                                 parameters => "w_crossover=bad\nw_width=bad\n");
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(G => [1, -1], F => [0.5, -0.25],
                      I => [0.2, [-0.3, -0.0675, 0.2, -0.8]]);
    unlink(map { "c-re$_.dat" } qw(G F I)) == 3 or die $!;

    is(run_sigmatrick(), 0, "new-im needs no real correlator files") or return;
    assert_result([map { 0.4 - 0.3 + $_ / 8 } @omega],
                  [-0.2375, -0.01, -0.01, -0.7375], "new-im");
    is_deeply([read_table("kk.input")], [read_table("imsigma.dat")],
              "KK receives the clipped actual ImSigma with no extra pi");
    assert_no_temporaries("successful KK");

    my @published = map { read_file($_) } @outputs;
    write_file("c-re$_.dat", "invalid unused correlator\n") for qw(G F I);
    is(run_sigmatrick(), 0, "new-im ignores malformed real correlators");
    is_deeply([map { read_file($_) } @outputs], \@published,
              "new-im outputs are independent of every real correlator");
};

subtest "crossover uses an even logarithmic quintic and clips only after blending" => sub {
    my $c = 0.5;
    my @mesh = map { $c * $_ } (-8, -4, -2, -1, -0.5, -0.25, -0.125,
                               0.125, 0.25, 0.5, 1, 2, 4, 8);
    my $dir = self_energy_fixture(
        clip => 0.01, sigma_h => 0.6, mode => "new-crossover", mesh => \@mesh,
        parameters => "w_crossover=$c\nw_width=4\n",
    );
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    # At r/c=0.5, 1, 2, t=1/4, 1/2, 3/4; these distinguish quintic from cubic.
    my %weight = (0.125 => 0, 0.25 => 0, 0.5 => 0.103515625,
                  1 => 0.5, 2 => 0.896484375, 4 => 1, 8 => 1);
    my @expected;
    for my $imi (-0.3, -0.05) {
        write_correlators(G => [1, -1], F => [0.5, -0.25], I => [0.2, $imi]);
        unlink("c-reI.dat") or die $!;
        @expected = map {
            my $w = $weight{abs($_) / $c};
            my $raw = (1 - $w) * ($imi + 0.0625) + $w * ($imi + 0.03125);
            $raw > -0.01 ? -0.01 : $raw;
        } @omega;
        is(run_sigmatrick(), 0, "crossover succeeds without reI (ImI=$imi)")
            or return;
        assert_result([map { 0.6 - 0.3 + $_ / 8 } @omega], \@expected,
                      "crossover ImI=$imi");
        is_deeply([read_table("kk.input")], [read_table("imsigma.dat")],
                  "crossover KK input agrees with the published clipped blend");
        assert_no_temporaries("crossover KK");
    }

    write_table("c-imG.dat", \@omega,
                [map { abs($_) >= $c * 4 ? 0 : 1 / $pi } @omega]);
    @expected = map {
        abs($omega[$_]) >= $c * 4 ? -0.01 : $expected[$_]
    } 0 .. $#omega;
    is(run_sigmatrick(), 0,
       "pure-new points, including the upper boundary, allow ImG=0 with G!=0");
    assert_result([map { 0.6 - 0.3 + $_ / 8 } @omega], \@expected,
                  "zero ImG in pure-new region");
};

subtest "imaginary-only division rejects zero ImG, including zero over zero" => sub {
    for my $case (["new-im", ""],
                  ["new-crossover", "w_crossover=1\nw_width=2\n"],
                  ["new-crossover", "w_crossover=0.25\nw_width=2\n"]) {
        my ($mode, $parameters) = @$case;
        my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.1,
                                     mode => $mode, parameters => $parameters);
        chdir($dir) or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        for my $imf (-0.25, 0) {
            write_correlators(G => [1, [-1, 0, -1, -1]],
                              F => [0.5, $imf], I => [0.2, -0.3]);
            write_output_sentinels();
            isnt(run_sigmatrick(), 0, "$mode rejects active ImG=0 with ImF=$imf");
            assert_output_sentinels("$mode singular input");
            assert_no_temporaries("$mode singular input");
        }
    }
};

subtest "mode and crossover parameters are strict" => sub {
    my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(G => [1, -1], F => [0.5, -0.25], I => [0.2, -0.3]);
    write_output_sentinels();
    my @invalid = (
        ["unknown mode", "sigmatrick=unknown\n"],
        ["empty mode", "sigmatrick=\n"],
        ["duplicate mode", "sigmatrick=new\nsigmatrick=new\n"],
        ["conflicting mode", "sigmatrick=new\nsigmatrick=old\n"],
    );
    for my $pair ([undef, 2], [1, undef], ["", 2], [0, 2], [-1, 2],
                  ["NaN", 2], ["1e999", 2], [1, ""], [1, 0], [1, 1],
                  [1, 0.9], [1, "Inf"], [1, "1e999"]) {
        my ($c, $b) = @$pair;
        my $parameters = "sigmatrick=new-crossover\n";
        $parameters .= "w_crossover=$c\n" if defined($c);
        $parameters .= "w_width=$b\n" if defined($b);
        my $label = "c=" . ($c // "absent") . ", b=" . ($b // "absent");
        push(@invalid, ["invalid crossover ($label)", $parameters]);
    }
    for my $case (@invalid) {
        write_file("param.loop", "clipSigma=0.01\n$case->[1]");
        isnt(run_sigmatrick(), 0, "$case->[0] is rejected");
        assert_output_sentinels($case->[0]);
    }
    ok(!-e "kk.args", "invalid parameters never reach KK");
};

subtest "KK modes require precisely their estimator inputs" => sub {
    for my $mode (qw(new-im new-crossover)) {
        my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.1, mode => $mode,
                                     parameters => "w_crossover=0.5\nw_width=2\n");
        chdir($dir) or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        write_correlators(G => [1, -1], F => [0.5, -0.25], I => [0.2, -0.3]);
        write_output_sentinels();
        for my $name (qw(imG imF imI), $mode eq "new-crossover" ? qw(reG reF) : ()) {
            rename("c-$name.dat", "saved.dat") or die $!;
            isnt(run_sigmatrick(), 0, "$mode requires $name");
            assert_output_sentinels("missing $name");
            rename("saved.dat", "c-$name.dat") or die $!;
        }
    }
};

subtest "KK failures preserve outputs and remove temporary files" => sub {
    for my $mode (qw(new-im new-crossover)) {
        my $dir = self_energy_fixture(clip => 0.01, sigma_h => 0.1, mode => $mode,
                                     parameters => "w_crossover=0.5\nw_width=2\n");
        chdir($dir) or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        write_correlators(G => [1, -1], F => [0.5, -0.25], I => [0.2, -0.3]);
        write_output_sentinels();
        for my $failure (qw(failure malformed nonfinite shifted truncated empty)) {
            local $ENV{SIGMA_KK_MODE} = $failure;
            unlink("kk.args") if -e "kk.args";
            isnt(run_sigmatrick(), 0, "$mode rejects $failure KK output");
            ok(-e "kk.args", "$failure reached the KK stub");
            assert_output_sentinels("$mode $failure");
            assert_no_temporaries("$mode $failure");
        }
    }
};

subtest "KK modes reject grids unsupported by the real backend" => sub {
    for my $mode (qw(new-im new-crossover)) {
        for my $mesh ([-1, 1], [-1, -0.25, 0, 0.25, 1], [-1, -0.2, 0.25, 1]) {
            my $dir = self_energy_fixture(
                clip => 0.01, sigma_h => 0.1, mode => $mode, mesh => $mesh,
                parameters => "w_crossover=0.5\nw_width=2\n",
            );
            chdir($dir) or die $!;
            local $ENV{PATH} = "$dir/bin:$original_path";
            write_correlators(G => [1, -1], F => [0.5, -0.25], I => [0.2, -0.3]);
            write_output_sentinels();
            isnt(run_sigmatrick(), 0, "$mode rejects grid @$mesh");
            assert_output_sentinels("unsupported KK grid");
            assert_no_temporaries("unsupported KK grid");
        }
    }
};

subtest "installed KK has the independent finite-support sign and normalization" => sub {
    if ($external_mode eq "0") {
        plan skip_all => "external numerical tools disabled by DMFT_TEST_EXTERNAL=0";
        return;
    }
    unless (grep { -x "$_/kk" } split(/:/, $original_path)) {
        die "missing external numerical tools: kk\n" if $external_mode eq "1";
        plan skip_all => "missing external numerical tools: kk";
        return;
    }
    for my $mode (qw(new-im new-crossover)) {
        my $dir = self_energy_fixture(
            clip => 0.01, sigma_h => 0.7, mode => $mode,
            mesh => [-2, -1, -0.25, 0.25, 1, 2],
            parameters => "w_crossover=0.5\nw_width=2\n",
        );
        chdir($dir) or die $!;
        unlink("bin/kk") or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        write_correlators(G => [1, -1], F => [0, 0], I => [9, -0.4]);
        unlink("c-reI.dat") or die $!;
        unlink("c-reG.dat", "c-reF.dat") == 2 or die $! if $mode eq "new-im";
        my @real = map {
            abs($_) == 2 ? 0.7 : 0.7 - 0.4 / $pi * log((2 - $_) / (2 + $_));
        } @omega;
        is(run_sigmatrick(), 0, "$mode uses the installed analytic Steffen KK")
            or next;
        assert_result(\@real, [(-0.4) x @omega], "$mode analytic KK");
        assert_no_temporaries("installed KK");
    }
};

subtest "equal-length shifted meshes are rejected before publication" => sub {
    my $dir = self_energy_fixture(clip => 1e-6, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(
        G => [1, -1],
        F => [0.5, -0.25],
        I => [0.2, -0.3],
    );
    write_table("c-imF.dat", [-1, -0.2, 0.25, 1], [(0.25 / $pi) x 4]);
    write_output_sentinels();

    isnt(run_sigmatrick(), 0, "shifted correlator mesh is rejected");
    assert_output_sentinels("mesh failure");

    write_correlators(
        G => [1, -1],
        F => [0.5, -0.25],
        I => [0.2, -0.3],
    );
    write_file("c-reI.dat", "-1 0\n-0.25 0\n0.25 1e999\n1 0\n");
    isnt(run_sigmatrick(), 0, "non-finite correlator value is rejected");
    assert_output_sentinels("non-finite input failure");
};

subtest "zero Green function is rejected before publication" => sub {
    my $dir = self_energy_fixture(clip => 1e-6, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(
        G => [0, 0],
        F => [0.5, -0.25],
        I => [0.2, -0.3],
    );
    write_output_sentinels();

    isnt(run_sigmatrick(), 0, "singular estimator input is rejected");
    assert_output_sentinels("singular input failure");
};

subtest "output preflight preserves an existing result set" => sub {
    my $dir = self_energy_fixture(clip => 1e-6, sigma_h => 0.1);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_correlators(
        G => [1, -1],
        F => [0.5, -0.25],
        I => [0.2, -0.3],
    );
    write_output_sentinels();
    unlink("imsigma.dat") or die $!;
    symlink("missing/imsigma.dat", "imsigma.dat") or die $!;

    isnt(run_sigmatrick(), 0, "invalid output target is rejected before writing");
    is(read_file("c-self.dat"), "preserve self\n",
       "publication failure preserves c-self.dat");
    is(read_file("resigma.dat"), "preserve real\n",
       "publication failure preserves resigma.dat");
    is(readlink("imsigma.dat"), "missing/imsigma.dat",
       "publication failure preserves the output symlink");

    unlink("imsigma.dat") or die $!;
    mkdir("imsigma.dat") or die $!;
    isnt(run_sigmatrick(), 0, "directory output target is rejected before writing");
    is(read_file("c-self.dat"), "preserve self\n",
       "rename preflight preserves c-self.dat");
    is(read_file("resigma.dat"), "preserve real\n",
       "rename preflight preserves resigma.dat");
    ok(-d "imsigma.dat", "rename preflight preserves the directory target");
};

chdir($original) or die $!;
done_testing();

sub self_energy_fixture {
    my %options = @_;
    my $dir = tempdir(CLEANUP => 1);
    make_path("$dir/bin");
    @omega = @{$options{mesh} // [-1, -0.25, 0.25, 1]};
    @gamma = @{$options{gamma} // [map {
        $_ == 0 || $_ == $#omega ? 0 : 0.1
    } 0 .. $#omega]};
    write_file("$dir/param.loop", "clipSigma=$options{clip}\n" .
               (exists($options{mode}) ? "sigmatrick=$options{mode}\n" : "") .
               ($options{parameters} // ""));
    write_file("$dir/param.mu.used", "0.2\n");
    write_file("$dir/customfdm.avg", "$options{sigma_h}\n");
    write_mesh("$dir/mesh.used.dat", \@omega);
    write_table("$dir/Delta.used.dat", \@omega, \@gamma);
    write_table("$dir/ReDelta.used.dat", \@omega, [(0.05) x @omega]);

    write_file("$dir/bin/extractcolumn", <<'EXTRACTCOLUMN');
#!/usr/bin/env perl
use strict;
use warnings;
my ($file, $column) = @ARGV;
$column eq "SigmaHd" or exit 1;
open(my $fh, "<", $file) or exit 1;
my $value = <$fh>;
defined($value) or exit 1;
print $value;
EXTRACTCOLUMN
    write_file("$dir/bin/kk", <<'KK');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy qw(copy);
@ARGV == 7 or die "KK requires explicit options and two temporary paths\n";
my ($input, $output) = splice(@ARGV, -2);
my %expected = ("-v" => undef, "--interpolation" => "steffen",
                "--algorithm" => "analytic");
while (@ARGV) {
    my $option = shift @ARGV;
    exists($expected{$option}) or die "Unexpected KK option: $option\n";
    my $value = delete($expected{$option});
    !defined($value) || shift(@ARGV) eq $value or die "Wrong KK option value\n";
}
!keys(%expected) or die "Missing KK options\n";
$input ne $output && $input !~ m{(?:^|/)imsigma\.dat\z}
    && $output !~ m{(?:^|/)resigma\.dat\z}
    or die "KK must use temporary input and output files\n";
open(my $log, ">", "kk.args") or die $!;
print {$log} "$input\n$output\n";
close($log) or die $!;
copy($input, "kk.input") or die $!;
open(my $in, "<", $input) or die $!;
my @rows = map { [split] } <$in>;
close($in) or die $!;
my $mode = $ENV{SIGMA_KK_MODE} // "success";
open(my $out, ">", $output) or die $!;
if ($mode eq "failure" || $mode eq "malformed") {
    print {$out} "malformed partial output\n";
} elsif ($mode ne "empty") {
    pop(@rows) if $mode eq "truncated";
    for my $row (@rows) {
        my $x = $row->[0] + ($mode eq "shifted" ? 0.01 : 0);
        my $value = $mode eq "nonfinite" ? "1e999" : $row->[0] / 8 - 0.3;
        print {$out} "$x $value\n";
    }
}
close($out) or die $!;
exit($mode eq "failure" ? 7 : 0);
KK
    chmod(0755, "$dir/bin/extractcolumn", "$dir/bin/kk");
    return $dir;
}

sub write_correlators {
    my %correlators = @_;
    for my $name (sort keys %correlators) {
        for my $part (0, 1) {
            my $value = $correlators{$name}[$part];
            my @values = ref($value) eq "ARRAY" ? @$value : ($value) x @omega;
            my $component = $part ? "im" : "re";
            write_table("c-$component$name.dat", \@omega,
                        [map { -$_ / $pi } @values]);
        }
    }
}

sub stage_production_layout {
    make_path("res");
    for my $file (qw(
        mesh.used.dat param.mu.used customfdm.avg Delta.used.dat ReDelta.used.dat
        c-reG.dat c-imG.dat c-reF.dat c-imF.dat c-reI.dat c-imI.dat
    )) {
        rename($file, "res/$file") or die "Can't stage $file: $!";
    }
    symlink("ReDelta.used.dat", "res/ReDelta.dat") or die $!;
}

sub run_sigmatrick {
    my ($result_dir) = @_;
    $result_dir //= ".";
    my $prefix = $result_dir eq "." ? "" : "$result_dir/";
    return system(
        $^X, $sigmatrick, $result_dir, "${prefix}mesh.used.dat",
        "${prefix}param.mu.used", "${prefix}customfdm.avg",
    );
}

sub write_output_sentinels {
    write_file("c-self.dat", "preserve self\n");
    write_file("imsigma.dat", "preserve imaginary\n");
    write_file("resigma.dat", "preserve real\n");
}

sub assert_output_sentinels {
    my ($context) = @_;
    is(read_file("c-self.dat"), "preserve self\n", "$context preserves c-self.dat");
    is(read_file("imsigma.dat"), "preserve imaginary\n",
       "$context preserves imsigma.dat");
    is(read_file("resigma.dat"), "preserve real\n",
        "$context preserves resigma.dat");
}

sub assert_result {
    my ($real, $imaginary, $context, $prefix) = @_;
    $prefix //= "";
    my @spectral = map {
        my $re = $omega[$_] + 0.2 - 0.05 - $real->[$_];
        my $im = $gamma[$_] - $imaginary->[$_];
        $im / ($pi * ($re**2 + $im**2));
    } 0 .. $#omega;
    for my $table (["resigma.dat", $real], ["imsigma.dat", $imaginary],
                   ["c-self.dat", \@spectral]) {
        my ($file, $expected) = @$table;
        my @actual = read_table("$prefix$file");
        is_deeply([map { $_->[0] } @actual], \@omega, "$context $file mesh");
        next if @actual != @omega;
        cmp_ok(max(map { abs($actual[$_][1] - $expected->[$_]) } 0 .. $#omega),
               "<=", 2e-13, "$context $file values");
    }
}

sub assert_no_temporaries {
    my ($context) = @_;
    my @temporary = (glob("*.tmp.*"), glob(".*.tmp.*"),
                     glob("res/*.tmp.*"), glob("res/.*.tmp.*"));
    if (-e "kk.args") {
        push(@temporary, grep { -e $_ || -l $_ } split(/\n/, read_file("kk.args")));
    }
    is_deeply(\@temporary, [], "$context removes temporary files");
}

sub write_table {
    my ($path, $mesh, $values) = @_;
    @$mesh == @$values or die "mesh/value length mismatch";
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    for my $index (0 .. $#$mesh) {
        printf {$fh} "%.17g %.17g\n", $mesh->[$index], $values->[$index];
    }
    close($fh) or die "Can't close $path: $!";
}

sub write_mesh {
    my ($path, $mesh) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print {$fh} "$_ ignored\n" for @$mesh;
    close($fh) or die "Can't close $path: $!";
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
        next if /^\s*(?:#.*)?$/;
        my @fields = split;
        @fields == 2 or die "Expected two columns in $path";
        push(@rows, [map { 0.0 + $_ } @fields]);
    }
    close($fh) or die "Can't close $path: $!";
    return @rows;
}
