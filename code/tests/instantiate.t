#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $scripts = "$FindBin::Bin/../scripts";
my $instantiate = "$scripts/instantiate";
my $wilson = "$scripts/wilson";
my $original = getcwd();
my $original_path = $ENV{PATH};

subtest "explicit Nmax remains supported" => sub {
    my $result = run_instantiate("Nmax=7\n");
    is($result->{status}, 0, "instantiate succeeds");
    assert_resolution($result, 7);
    is($result->{copied_param}, $result->{source_param},
       "the solver receives the original parameter file");
};

subtest "Tmin follows nrginit for every checked-in twist" => sub {
    my @cases = (
        [0.25, 54],
        [0.50, 54],
        [0.75, 53],
        [1.00, 53],
    );
    for my $case (@cases) {
        my ($z, $expected) = @$case;
        my $result = run_instantiate(join("\n",
            "Tmin=1e-7",
            "Lambda=2",
            "z=$z",
            "bandrescale=10",
            "discretization=Z",
            "",
        ));
        is($result->{status}, 0, "z=$z succeeds");
        assert_resolution($result, $expected, "z=$z");
        unlike($result->{copied_param} // "", qr/^Nmax=/m,
               "z=$z copied param does not expose synthetic Nmax");
    }
};

subtest "the nrginit threshold comparison is inclusive" => sub {
    my $result = run_instantiate(join("\n",
        "Tmin=0.15625",
        "Lambda=4",
        "z=1",
        "bandrescale=1",
        "discretization=Y",
    ));
    is($result->{status}, 0, "exact-threshold calculation succeeds");
    assert_resolution($result, 3);

    my $rounding = run_instantiate(join("\n",
        "Tmin=0.34756349141427556",
        "Lambda=1.1",
        "z=0.1",
        "bandrescale=1",
        "discretization=Y",
        "",
    ));
    is($rounding->{status}, 0, "nrginit evaluation-order case succeeds");
    assert_resolution($rounding, 24, "nrginit evaluation order");
};

subtest "Tmin_ratio and explicit Tmin use nrginit precedence" => sub {
    my $ratio = run_instantiate(join("\n",
        "Nmax=999",
        "T=0.3125",
        "Tmin_ratio=0.5",
        "Lambda=4",
        "z=1",
        "bandrescale=1",
        "discretization=Y",
        "",
    ));
    is($ratio->{status}, 0, "positive ratio calculation succeeds");
    assert_resolution($ratio, 3, "positive ratio");
    like($ratio->{copied_param} // "", qr/^Nmax=999$/m,
         "synthetic resolution does not rewrite the copied Nmax");

    my $explicit = run_instantiate(join("\n",
        "T=1",
        "Tmin_ratio=1",
        "Tmin=0.15625",
        "Lambda=4",
        "z=1",
        "bandrescale=1",
        "discretization=Y",
        "",
    ));
    is($explicit->{status}, 0, "explicit Tmin calculation succeeds");
    assert_resolution($explicit, 3, "explicit Tmin");

    my $inactive = run_instantiate("Nmax=7\nT=1\nTmin_ratio=0\n");
    is($inactive->{status}, 0, "non-positive ratio leaves Nmax active");
    assert_resolution($inactive, 7, "non-positive ratio");
};

subtest "invalid cutoff configurations fail before nrgchain" => sub {
    my @cases = (
        ["Nmax=7\nTmin=1e-3\n", "explicit Nmax and Tmin"],
        ["Tmin=invalid\n", "invalid Tmin"],
        ["Nmax=7\nT=1e-9999\nTmin_ratio=1\n",
         "underflowing ratio temperature"],
        ["Nmax=999\n", "Nmax outside nrginit's range"],
        ["Tmin=1e-3\nLambda=1\n", "invalid Lambda"],
        ["Tmin=2\nLambda=4\nz=1\nbandrescale=1\ndiscretization=Y\n",
         "cutoff above first Wilson scale"],
        ["Tmin_ratio=1\n", "ratio without T or Nmax"],
    );
    for my $case (@cases) {
        my ($settings, $name) = @$case;
        my $result = run_instantiate($settings);
        isnt($result->{status}, 0, "$name is rejected");
        ok(!defined($result->{chain_param}), "$name fails before nrgchain");
    }
};

chdir($original) or die $!;
done_testing();

sub run_instantiate {
    my $settings = shift;
    my $dir = tempdir(CLEANUP => 1);
    make_path("$dir/bin", "$dir/scripts", "$dir/template", "$dir/out");
    symlink($wilson, "$dir/scripts/wilson") or die "Can't link wilson: $!";

    my $source_param = "[extra]\n[dmft]\n[param]\npolarized=false\n" . $settings;
    write_file("$dir/param", $source_param);
    write_file("$dir/Delta.dat", "-1 1\n1 1\n");
    write_file("$dir/template/data.in", "1 0 0\n");
    write_file("$dir/bin/nrgchain", <<'NRGCHAIN');
#!/usr/bin/env perl
use strict;
use warnings;
my $parameter_file = shift // "param";
@ARGV and die "unexpected nrgchain arguments\n";
open(my $input, "<", $parameter_file) or die $!;
local $/;
my $parameters = <$input>;
close($input) or die $!;
my ($nmax) = $parameters =~ /^\s*Nmax\s*=\s*(\d+)\s*$/m;
defined($nmax) or die "Nmax missing from nrgchain input\n";
open(my $capture, ">", "nrgchain.input") or die $!;
print {$capture} $parameters;
close($capture) or die $!;
for my $filename (qw(xi.dat zeta.dat)) {
    open(my $output, ">", $filename) or die $!;
    print {$output} "$_\n" for 0 .. $nmax;
    close($output) or die $!;
}
open(my $theta, ">", "theta.dat") or die $!;
print {$theta} "1\n";
close($theta) or die $!;
NRGCHAIN
    chmod(0755, "$dir/bin/nrgchain") or die "Can't chmod nrgchain: $!";

    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    my ($writer, $reader);
    my $error = gensym;
    my $pid = open3($writer, $reader, $error, $^X, $instantiate, "out");
    close($writer);
    local $/;
    my $stdout = <$reader> // "";
    my $stderr = <$error> // "";
    waitpid($pid, 0);
    my $status = $?;

    my %result = (
        status => $status,
        stdout => $stdout,
        stderr => $stderr,
        source_param => $source_param,
    );
    $result{chain_param} = read_file("nrgchain.input") if -e "nrgchain.input";
    $result{data} = read_file("out/data") if -e "out/data";
    $result{copied_param} = read_file("out/param") if -e "out/param";
    chdir($original) or die $!;
    return \%result;
}

sub assert_resolution {
    my ($result, $expected, $context) = @_;
    $context //= "resolution";
    ok(defined($result->{chain_param}), "$context reaches nrgchain");
    return if !defined($result->{chain_param});

    my @values = $result->{chain_param} =~ /^\s*Nmax\s*=\s*(\d+)\s*$/mg;
    is(scalar(@values), 1, "$context passes exactly one Nmax to nrgchain");
    is($values[0], $expected, "$context passes the expected Nmax") if @values;
    like($result->{data} // "", qr/^1\s+$expected\s+0$/m,
         "$context writes the same Nmax to data");
    my @data_lines = split(/\n/, $result->{data} // "");
    is(scalar(@data_lines), 2*$expected + 6,
       "$context writes Nmax+1 coefficients for both tables");
    is($data_lines[2], $expected, "$context labels the xi table");
    is($data_lines[$expected + 4], $expected,
       "$context labels the zeta table after Nmax+1 xi values");
}

sub write_file {
    my ($path, $contents) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print {$fh} $contents;
    close($fh) or die "Can't close $path: $!";
}

sub read_file {
    my $path = shift;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    local $/;
    my $contents = <$fh>;
    close($fh) or die "Can't close $path: $!";
    return $contents;
}
