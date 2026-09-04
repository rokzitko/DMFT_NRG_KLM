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

my $validator = "$FindBin::Bin/../scripts/validate_bare_inputs";
my $asym = "$FindBin::Bin/../../more_examples/asym";
my $original = getcwd();
my $original_path = $ENV{PATH};

sub write_file {
    my ($path, $contents) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print {$fh} $contents;
    close($fh) or die "Can't close $path: $!";
}

sub fixture {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("bin");
    write_file("DOS.dat", "-1 0\n0 1\n1 0\n");
    write_file("PHI.dat", "-1 0\n0 0.75\n1 0\n");
    write_file("param.eps", "0.25\n");
    write_file(
        "bin/integ",
        <<'STUB',
#!/usr/bin/env perl
use strict;
use warnings;
open(my $log, ">>", $ENV{INTEG_LOG}) or die $!;
print {$log} join(" ", @ARGV), "\n";
close($log) or die $!;
if ($ENV{INTEG_FAIL}) {
    print STDERR "stub integration failure\n";
    exit 7;
}
my $file = $ARGV[-1];
if (grep { $_ eq "--energy-moment" } @ARGV) {
    print(($ENV{DOS_MOMENT} // "0.25"), "\n");
} elsif ($file =~ /PHI/) {
    print(($ENV{PHI_TOTAL} // "0.75"), "\n");
} else {
    print(($ENV{DOS_TOTAL} // "1"), "\n");
}
STUB
    );
    chmod(0755, "bin/integ");
    $ENV{INTEG_LOG} = "$dir/integ.log";
    $ENV{PATH} = "$dir/bin:$original_path";
    delete @ENV{qw(INTEG_FAIL DOS_TOTAL DOS_MOMENT PHI_TOTAL)};
    return $dir;
}

sub run_validator {
    my (@arguments) = @_;
    my $error = gensym();
    my $pid = open3(undef, my $output, $error, $^X, $validator, @arguments);
    local $/;
    my $stdout = <$output> // "";
    my $stderr = <$error> // "";
    waitpid($pid, 0);
    return ($? >> 8, $stdout, $stderr);
}

subtest "valid represented tables" => sub {
    my $dir = fixture();
    my ($status, $stdout, $stderr) = run_validator("--require-phi");
    is($status, 0, "valid DOS, centroid, and transport mesh are accepted");
    like($stdout, qr/DOS integral=1, centroid=0\.25, param\.eps=0\.25/,
         "validation summary reports the represented moments");
    like($stdout, qr/PHI integral=0\.75/, "summary reports the transport weight");
    is($stderr, "", "valid input emits no diagnostics");

    my @commands = split(/\n/, read_file("$dir/integ.log"));
    is(scalar(@commands), 3, "DOS total, DOS moment, and PHI total are integrated");
    like($commands[0], qr/^--interpolation steffen --total -- DOS\.dat$/,
         "DOS weight uses the Steffen representation");
    like($commands[1], qr/^--interpolation steffen --energy-moment -- DOS\.dat$/,
         "DOS first moment uses the Steffen representation");
    like($commands[2], qr/^--interpolation steffen --total -- PHI\.dat$/,
         "transport weight uses the Steffen representation");
};

subtest "normalization and centroid tolerances" => sub {
    fixture();
    local $ENV{DOS_TOTAL} = "1.000000005";
    local $ENV{DOS_MOMENT} = "0.25000000125";
    my ($status) = run_validator();
    is($status, 0, "roundoff-scale normalization error is accepted");

    local $ENV{DOS_TOTAL} = "1.00000002";
    my ($bad_status, undef, $bad_stderr) = run_validator();
    isnt($bad_status, 0, "material normalization error is rejected");
    like($bad_stderr, qr/DOS normalization is not one/,
         "normalization failure is explicit");

    local $ENV{DOS_TOTAL} = "1";
    write_file("param.eps", "0.2500000002\n");
    my ($moment_status, undef, $moment_stderr) = run_validator();
    isnt($moment_status, 0, "inconsistent static shift is rejected");
    like($moment_stderr, qr/param\.eps does not match.*first moment/,
         "centroid failure is explicit");
};

subtest "strict table validation" => sub {
    fixture();
    write_file("DOS.dat", "# epsilon rho\n-1 0\n0 1\n1 0\n");
    write_file("PHI.dat", "  # epsilon Phi\n-1 0\n0 0.75\n1 0\n");
    my ($comment_status) = run_validator();
    is($comment_status, 0, "whole-line comments accepted by integ are accepted");

    write_file("DOS.dat", "-1 0\n0 -0.1\n1 0\n");
    my ($negative_status, undef, $negative_stderr) = run_validator();
    isnt($negative_status, 0, "negative DOS is rejected before integration");
    like($negative_stderr, qr/negative DOS value.*line 2/,
         "negative row is identified");

    write_file("DOS.dat", "-1 0\n0 1\n0 0.5\n1 0\n");
    my ($mesh_status, undef, $mesh_stderr) = run_validator();
    isnt($mesh_status, 0, "duplicate DOS knot is rejected");
    like($mesh_stderr, qr/mesh is not strictly increasing.*line 3/,
         "bad knot is identified");

    write_file("DOS.dat", "-1 0\n0 1e-9999\n1 0\n");
    my ($underflow_status, undef, $underflow_stderr) = run_validator();
    isnt($underflow_status, 0, "underflowed DOS value is rejected");
    like($underflow_stderr, qr/numeric underflow.*line 2/,
         "underflow is distinguished from an exact zero");

    write_file("DOS.dat", "-1 0\n0 1\n1 0\n");
    write_file("param.eps", "nan\n");
    my ($scalar_status, undef, $scalar_stderr) = run_validator();
    isnt($scalar_status, 0, "non-finite param.eps is rejected");
    like($scalar_stderr, qr/invalid numeric value.*param\.eps/,
         "invalid scalar is identified");
};

subtest "asymmetric generators use the validation contract" => sub {
    fixture();
    unlink("PHI.dat") or die $!;
    local $ENV{DOS_MOMENT} = "0.175";
    is(system("$asym/mkDOS"), 0, "asymmetric generator creates a DOS");
    is(system("$asym/mkparameps"), 0,
       "asymmetric moment helper creates param.eps");
    my ($status) = run_validator();
    is($status, 0, "generated files pass the common validator");
    my @rows = grep { /\S/ } split(/\n/, read_file("DOS.dat"));
    is(scalar(@rows), 2601, "generator emits the edge-refined table size");
    is(0.0 + read_file("param.eps"), 0.175,
       "moment helper publishes the represented centroid");
};

subtest "optional and matching transport table" => sub {
    fixture();
    unlink("PHI.dat") or die $!;
    my ($optional_status, $optional_stdout) = run_validator();
    is($optional_status, 0, "core DMFT accepts an absent PHI table");
    like($optional_stdout, qr/PHI absent/, "absence is reported");

    my ($required_status, undef, $required_stderr) =
        run_validator("--require-phi");
    isnt($required_status, 0, "transport mode requires PHI.dat");
    like($required_stderr, qr/required transport table PHI\.dat is missing/,
         "missing transport input is explicit");

    write_file("PHI.dat", "-1 0\n0.1 0.75\n1 0\n");
    my ($mesh_status, undef, $mesh_stderr) = run_validator();
    isnt($mesh_status, 0, "a stale PHI mesh is rejected");
    like($mesh_stderr, qr/frequency mesh mismatch.*row 2/,
         "transport mesh mismatch is identified");

    write_file("PHI.dat", "-1 0\n0 0.75\n1 0\n");
    local $ENV{PHI_TOTAL} = "0";
    my ($weight_status, undef, $weight_stderr) = run_validator();
    isnt($weight_status, 0, "zero represented transport weight is rejected");
    like($weight_stderr, qr/PHI integral must be positive/,
         "transport weight failure is explicit");
};

subtest "integration backend failures are fatal" => sub {
    fixture();
    local $ENV{INTEG_FAIL} = 1;
    my ($status, undef, $stderr) = run_validator();
    isnt($status, 0, "integ failure rejects the inputs");
    like($stderr, qr/stub integration failure.*integ failed/s,
         "backend diagnostic remains visible");
};

chdir($original) or die $!;
done_testing();

sub read_file {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    local $/;
    my $contents = <$fh>;
    close($fh) or die "Can't close $path: $!";
    return $contents;
}
