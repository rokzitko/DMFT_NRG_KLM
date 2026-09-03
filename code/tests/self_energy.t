#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use Test::More;

my $sigmatrick = "$FindBin::Bin/../scripts/sigmatrick";
my $original = getcwd();
my $original_path = $ENV{PATH};
my $pi = 4 * atan2(1, 1);
my @omega = (-1, -0.25, 0.25, 1);
my @gamma = (0, 0.1, 0.1, 0);

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
    write_file("$dir/param.loop", "clipSigma=$options{clip}\n");
    write_file("$dir/param.mu.used", "0.2\n");
    write_file("$dir/customfdm.avg", "$options{sigma_h}\n");
    write_mesh("$dir/mesh.used.dat", \@omega);
    write_table("$dir/Delta.used.dat", \@omega, \@gamma);
    write_table("$dir/ReDelta.used.dat", \@omega, [(0.05) x @omega]);

    write_file("$dir/bin/getparam", <<'GETPARAM');
#!/usr/bin/env perl
use strict;
use warnings;
my ($name, $file) = @ARGV;
open(my $fh, "<", $file) or exit 1;
while (<$fh>) {
    if (/^\s*\Q$name\E\s*=\s*(.*?)\s*$/) {
        print "$1\n";
        exit 0;
    }
}
exit 1;
GETPARAM
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
    chmod(0755, "$dir/bin/getparam", "$dir/bin/extractcolumn");
    return $dir;
}

sub write_correlators {
    my %correlators = @_;
    for my $name (qw(G F I)) {
        my ($real, $imaginary) = @{$correlators{$name}};
        write_table("c-re$name.dat", \@omega, [(-$real / $pi) x @omega]);
        write_table("c-im$name.dat", \@omega, [(-$imaginary / $pi) x @omega]);
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
