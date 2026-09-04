#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use FindBin;
use IPC::Open3 qw(open3);
use Symbol qw(gensym);
use Test::More;

my $causal_delta = "$FindBin::Bin/../scripts/causalDelta";
my $original = getcwd();
my $original_path = $ENV{PATH};
my $safe_gamma =
    "-2 1e-9\n-1 0.05\n-0.5 0.2\n0.5 -5e-8\n1 0.4\n2 -2e-9\n";

subtest "bounded numerical corrections produce a causal triplet" => sub {
    my $dir = causality_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file("Gamma.raw.dat", $safe_gamma);

    my $result = run_causal_delta();
    is($result->{status}, 0, "bounded projection succeeds");
    like($result->{stdout}, qr/raw_min=-5e-08\b/, "raw minimum is reported");
    like($result->{stdout}, qr/endpoint_policy=strict\b/,
         "strict endpoint policy is reported");
    like($result->{stdout}, qr/negative_points=2\b/,
         "negative correction count is reported");
    like($result->{stdout}, qr/endpoint_points=2\b/,
         "endpoint correction count is reported");
    like($result->{stdout}, qr/floored_points=2\b/,
         "interior floor count is reported");
    my ($correction_l1) = $result->{stdout} =~ /correction_l1=([0-9.eE+-]+)/;
    cmp_ok(abs($correction_l1 - 0.112500039), "<=", 2e-15,
           "integrated correction size is reported accurately");

    my @gamma = read_table("Delta.dat");
    my @imaginary = read_table("ImDelta.dat");
    my @real = read_table("ReDelta.dat");
    is_deeply([map { $_->[1] } @gamma], [0, 0.1, 0.2, 0.1, 0.4, 0],
              "small corrections are projected to the documented floor and guards");
    for my $index (0 .. $#gamma) {
        cmp_ok(abs($imaginary[$index][1] + $gamma[$index][1]), "<=", 1e-16,
               "ImDelta equals -Gamma at row " . ($index + 1));
        cmp_ok(abs($real[$index][1] - 0.25), "<=", 1e-16,
               "ReDelta contains the static shift at row " . ($index + 1));
    }
    my @temporary = glob("*.tmp.*");
    is(scalar(@temporary), 0, "successful projection leaves no temporary files");
};

subtest "materially negative Gamma is rejected without publication" => sub {
    my $dir = causality_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    publish_safe_triplet();
    my %published = published_outputs();
    write_file(
        "Gamma.raw.dat",
        "-2 0\n-1 0.05\n-0.5 0.2\n0.5 -1e-3\n1 0.4\n2 0\n",
    );

    my $result = run_causal_delta();
    isnt($result->{status}, 0, "material negative value is rejected");
    like($result->{stderr}, qr/below the allowed negative tolerance/,
         "negative-value failure is descriptive");
    assert_outputs_unchanged(\%published, "negative-value failure");
};

subtest "material endpoint tails are rejected without publication" => sub {
    my $dir = causality_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    publish_safe_triplet();
    my %published = published_outputs();
    write_file(
        "Gamma.raw.dat",
        "-2 1e-4\n-1 0.05\n-0.5 0.2\n0.5 0.2\n1 0.4\n2 0\n",
    );

    my $result = run_causal_delta();
    isnt($result->{status}, 0, "material endpoint is rejected");
    like($result->{stderr}, qr/endpoint Gamma .* exceeds tolerance/,
         "endpoint failure is descriptive");
    assert_outputs_unchanged(\%published, "endpoint failure");
};

subtest "legacy mode accepts only floor-sized endpoint guards" => sub {
    my $dir = causality_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file(
        "Gamma.raw.dat",
        "-2 0.1\n-1 0.2\n-0.5 0.2\n0.5 0.2\n1 0.2\n2 -0.1\n",
    );

    my $strict = run_causal_delta();
    isnt($strict->{status}, 0, "strict mode rejects floor-sized endpoints");
    my $legacy = run_causal_delta(allow_floor_endpoints => 1);
    is($legacy->{status}, 0, "legacy mode accepts floor-sized endpoints");
    like($legacy->{stdout}, qr/endpoint_policy=floor\b/,
         "legacy endpoint policy is reported");
    my @gamma = read_table("Delta.dat");
    is($gamma[0][1], 0, "legacy lower endpoint is canonicalized");
    is($gamma[-1][1], 0, "legacy upper endpoint is canonicalized");

    write_file(
        "Gamma.raw.dat",
        "-2 0.1001\n-1 0.2\n-0.5 0.2\n0.5 0.2\n1 0.2\n2 0\n",
    );
    my $oversized = run_causal_delta(allow_floor_endpoints => 1);
    isnt($oversized->{status}, 0, "legacy mode rejects endpoints above the floor");
};

subtest "peak-relative endpoint tolerance is enforced at production scale" => sub {
    my $dir = causality_fixture(clip => 1e-6);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file(
        "Gamma.raw.dat",
        "-2 -3e-9\n-1 0.1\n-0.5 0.2\n0.5 0.2\n1 0.4\n2 3e-9\n",
    );
    is(run_causal_delta()->{status}, 0,
       "endpoint below peak-relative tolerance is accepted for either sign");

    write_file(
        "Gamma.raw.dat",
        "-2 5e-9\n-1 0.1\n-0.5 0.2\n0.5 0.2\n1 0.4\n2 0\n",
    );
    isnt(run_causal_delta()->{status}, 0,
         "endpoint above peak-relative tolerance is rejected");
};

subtest "peak-relative negative tolerance is enforced" => sub {
    my $dir = causality_fixture(clip => 1e-15);
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file(
        "Gamma.raw.dat",
        "-2 0\n-1 0.1\n-0.5 0.2\n0.5 -3e-13\n1 0.4\n2 0\n",
    );
    is(run_causal_delta()->{status}, 0,
       "negative value below peak-relative tolerance is accepted");

    write_file(
        "Gamma.raw.dat",
        "-2 0\n-1 0.1\n-0.5 0.2\n0.5 -5e-13\n1 0.4\n2 0\n",
    );
    isnt(run_causal_delta()->{status}, 0,
         "negative value above peak-relative tolerance is rejected");
};

subtest "projection cannot manufacture a bath from nonpositive input" => sub {
    my $dir = causality_fixture();
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    publish_safe_triplet();
    my %published = published_outputs();
    write_file("Gamma.raw.dat", "-2 0\n-1 0\n-0.5 0\n0.5 0\n1 0\n2 0\n");

    my $result = run_causal_delta();
    isnt($result->{status}, 0, "nonpositive bath is rejected");
    like($result->{stderr}, qr/no positive interior Gamma/,
         "missing-positive-bath failure is descriptive");
    assert_outputs_unchanged(\%published, "nonpositive-bath failure");
};

chdir($original) or die $!;
done_testing();

sub causality_fixture {
    my %options = @_;
    my $clip = $options{clip} // 0.1;
    my $dir = tempdir(CLEANUP => 1);
    make_path("$dir/bin");
    write_file("$dir/param.loop", "clipDelta=$clip\n");
    write_file("$dir/param.eps", "0.25\n");
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
    write_file("$dir/bin/kk", <<'KK');
#!/usr/bin/env perl
use strict;
use warnings;
my ($input, $output) = @ARGV[-2, -1];
open(my $in, "<", $input) or die $!;
open(my $out, ">", $output) or die $!;
while (<$in>) {
    next if /^\s*(?:#.*)?$/;
    my ($omega) = split;
    print {$out} "$omega 0\n";
}
close($in) or die $!;
close($out) or die $!;
KK
    chmod(0755, "$dir/bin/getparam", "$dir/bin/kk");
    return $dir;
}

sub publish_safe_triplet {
    write_file("Gamma.raw.dat", $safe_gamma);
    my $result = run_causal_delta();
    $result->{status} == 0 or die "Unable to publish test fixture: $result->{stderr}";
}

sub run_causal_delta {
    my %options = @_;
    my $error = gensym;
    my ($writer, $reader);
    my @arguments;
    push(@arguments, "--allow-floor-endpoints")
        if $options{allow_floor_endpoints};
    my $pid = open3(
        $writer, $reader, $error, $^X, $causal_delta,
        @arguments,
        "Gamma.raw.dat", "Delta.dat", "ReDelta.dat", "ImDelta.dat",
    );
    close($writer);
    my ($stdout, $stderr);
    {
        local $/;
        $stdout = <$reader> // "";
        $stderr = <$error> // "";
    }
    close($reader);
    close($error);
    waitpid($pid, 0);
    return {status => $?, stdout => $stdout, stderr => $stderr};
}

sub published_outputs {
    return map { $_ => read_file($_) } qw(Delta.dat ReDelta.dat ImDelta.dat);
}

sub assert_outputs_unchanged {
    my ($published, $context) = @_;
    for my $file (qw(Delta.dat ReDelta.dat ImDelta.dat)) {
        is(read_file($file), $published->{$file}, "$context preserves $file");
    }
    my @temporary = glob("*.tmp.*");
    is(scalar(@temporary), 0, "$context leaves no temporary files");
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
