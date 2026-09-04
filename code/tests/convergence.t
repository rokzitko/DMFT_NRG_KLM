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

my $checkconv = "$FindBin::Bin/../scripts/checkconv";
my $original = getcwd();
my $original_path = $ENV{PATH};

subtest "minimum and maximum iterations are inclusive" => sub {
    my $dir = convergence_fixture(<<'PARAMETERS');
conveps=1e-4
convwindow=2
miniter=3
maxiter=5
PARAMETERS
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";

    write_file("DIFFS_C", "2 1e-6\n");
    is(run_checkconv(2), 0, "iteration below miniter produces a decision");
    is(read_file("DECISION"), "continue\n", "iteration below miniter continues");

    write_file("DIFFS_C", "2 1e-6\n3 2e-6\n");
    is(run_checkconv(3), 0, "iteration equal to miniter produces a decision");
    is(read_file("DECISION"), "converged\n", "convergence is allowed at miniter");

    write_file("DIFFS_C", "2 1e-2\n3 1e-2\n4 1e-6\n5 2e-6\n");
    is(run_checkconv(5), 0, "converged maxiter cycle produces a decision");
    is(read_file("DECISION"), "converged\n",
       "convergence takes precedence at maxiter");

    write_file("DIFFS_C", "2 1e-2\n3 1e-2\n4 1e-6\n5 2e-3\n");
    is(run_checkconv(5), 0, "nonconverged maxiter cycle produces a decision");
    is(read_file("DECISION"), "stop\n", "iteration equal to maxiter stops");

    write_file(
        "DIFFS_C",
        "2 1e-6\n3 1e-6\n4 1e-6\n5 1e-6\n6 1e-6\n",
    );
    is(run_checkconv(6), 0, "cycle beyond maxiter produces a decision");
    is(read_file("DECISION"), "stop\n", "cycle beyond maxiter cannot converge");
};

subtest "convergence uses the maximum over the complete window" => sub {
    my $dir = convergence_fixture(<<'PARAMETERS');
conveps=1e-4
convwindow=3
miniter=2
maxiter=10
PARAMETERS
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";

    write_file("DIFFS_C", "2 1e-6\n3 1e-2\n4 2e-6\n");
    is(run_checkconv(4), 0, "window with one large residual produces a decision");
    is(read_file("DECISION"), "continue\n",
       "one large residual blocks convergence");

    write_file("DIFFS_C", "2 1e-6\n3 3e-6\n4 2e-6\n");
    is(run_checkconv(4), 0, "closed window produces a decision");
    is(read_file("DECISION"), "converged\n",
       "all rows below tolerance converge");
};

subtest "an incomplete window stops at the iteration limit" => sub {
    my $dir = convergence_fixture(<<'PARAMETERS');
conveps=1e-4
convwindow=3
miniter=1
maxiter=3
PARAMETERS
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";

    write_file("DIFFS_C", "2 1e-6\n");
    is(run_checkconv(2), 0, "incomplete pre-limit window produces a decision");
    is(read_file("DECISION"), "continue\n", "incomplete pre-limit window continues");

    write_file("DIFFS_C", "2 1e-6\n3 1e-6\n");
    is(run_checkconv(3), 0, "incomplete maxiter window produces a decision");
    is(read_file("DECISION"), "stop\n", "incomplete window cannot exceed maxiter");
};

subtest "missing optional occupancy limits use silent defaults" => sub {
    my $dir = convergence_fixture(valid_parameters());
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file("DIFFS_C", "2 1e-3\n");

    my $result = run_checkconv_captured(2);
    is($result->{status}, 0, "missing optional limits produce a decision");
    is($result->{stderr}, "", "optional getparam diagnostics are suppressed");
    is(read_file("DECISION"), "continue\n", "default limits permit evaluation");
};

subtest "invalid convergence histories preserve the prior decision" => sub {
    my @cases = (
        ["2 1e-6\n4 1e-6\n", 4, "iteration gap"],
        ["2 1e-6\n3 1e-6\n", 4, "stale final iteration"],
        ["2 1e999\n", 2, "non-finite difference"],
        ["2 -1e-6\n", 2, "negative difference"],
        ["2 broken\n", 2, "nonnumeric difference"],
        ["2 1e-6 extra\n", 2, "extra column"],
        ["2 1e-6\n2 1e-6\n", 2, "duplicate iteration"],
        ["1 1e-6\n2 1e-6\n", 2, "iteration-one difference"],
        ["2 1e-6\n", 1, "history present at iteration one"],
    );

    for my $case (@cases) {
        my ($history, $iteration, $name) = @$case;
        my $dir = convergence_fixture(valid_parameters());
        chdir($dir) or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        write_file("DIFFS_C", $history);
        write_file("DECISION", "preserve\n");
        isnt(run_checkconv($iteration), 0, "$name is rejected");
        is(read_file("DECISION"), "preserve\n", "$name preserves DECISION");
    }

    my $dir = convergence_fixture(valid_parameters());
    chdir($dir) or die $!;
    local $ENV{PATH} = "$dir/bin:$original_path";
    write_file("DECISION", "preserve\n");
    isnt(run_checkconv(2), 0, "missing history is rejected after iteration one");
    is(read_file("DECISION"), "preserve\n", "missing history preserves DECISION");
};

subtest "invalid convergence parameters preserve the prior decision" => sub {
    my @cases = (
        ["conveps=0\nconvwindow=1\nminiter=1\nmaxiter=10\n", "zero conveps"],
        ["conveps=1e999\nconvwindow=1\nminiter=1\nmaxiter=10\n",
         "non-finite conveps"],
        ["conveps=1e-4\nconvwindow=1.5\nminiter=1\nmaxiter=10\n",
         "fractional convwindow"],
        ["conveps=1e-4\nconvwindow=0\nminiter=1\nmaxiter=10\n",
         "zero convwindow"],
        ["conveps=1e-4\nconvwindow=1\nminiter=-1\nmaxiter=10\n",
         "negative miniter"],
        ["conveps=1e-4\nconvwindow=1\nminiter=1\nmaxiter=0\n",
         "zero maxiter"],
        ["conveps=1e-4\nconvwindow=1\nminiter=4\nmaxiter=3\n",
         "maxiter below miniter"],
    );

    for my $case (@cases) {
        my ($parameters, $name) = @$case;
        my $dir = convergence_fixture($parameters);
        chdir($dir) or die $!;
        local $ENV{PATH} = "$dir/bin:$original_path";
        write_file("DIFFS_C", "2 1e-6\n");
        write_file("DECISION", "preserve\n");
        isnt(run_checkconv(2), 0, "$name is rejected");
        is(read_file("DECISION"), "preserve\n", "$name preserves DECISION");
    }
};

chdir($original) or die $!;
done_testing();

sub convergence_fixture {
    my ($parameters) = @_;
    my $dir = tempdir(CLEANUP => 1);
    make_path("$dir/bin");
    write_file("$dir/param.loop", $parameters);
    write_file("$dir/bin/getparam", <<'GETPARAM');
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
print STDERR "missing parameter $name\n";
exit 1;
GETPARAM
    chmod(0755, "$dir/bin/getparam");
    return $dir;
}

sub valid_parameters {
    return "conveps=1e-4\nconvwindow=1\nminiter=1\nmaxiter=10\n";
}

sub run_checkconv {
    my ($iteration) = @_;
    return system(
        $^X, $checkconv, "--iteration", $iteration,
        "--diffs", "DIFFS_C", "--decision", "DECISION",
    );
}

sub run_checkconv_captured {
    my ($iteration) = @_;
    my $error = gensym;
    my ($writer, $reader);
    my $pid = open3(
        $writer, $reader, $error,
        $^X, $checkconv, "--iteration", $iteration,
        "--diffs", "DIFFS_C", "--decision", "DECISION",
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
