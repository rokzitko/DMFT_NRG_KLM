#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd);
use File::Temp qw(tempdir);
use FindBin;
use Test::More;

my $fixture = "$FindBin::Bin/real_nrg";
my $original = getcwd();

sub write_file {
    my ($path, $contents) = @_;
    open(my $fh, ">", $path) or die "Can't write $path: $!";
    print {$fh} $contents;
    close($fh) or die "Can't close $path: $!";
}

subtest "seed is symmetric and covers the configured support" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    is(system("python3", "$fixture/make_seed.py", "Delta.next.dat", "mesh.next.dat"),
       0, "seed generator succeeds");
    my @gamma = read_table("Delta.next.dat");
    my @mesh = read_table("mesh.next.dat");
    is(scalar(@gamma), scalar(@mesh), "seed and mesh have the same size");
    cmp_ok(scalar(@gamma), ">", 500, "reduced broadening mesh remains resolved");
    is($gamma[0][0], -40, "lower support guard is exact");
    is($gamma[-1][0], 40, "upper support guard is exact");
    is($gamma[0][1], 0, "lower support value is zero");
    is($gamma[-1][1], 0, "upper support value is zero");
    ok(grep({ $_->[0] == -1 && $_->[1] == 0 } @gamma),
       "negative band edge is explicit");
    ok(grep({ $_->[0] == 1 && $_->[1] == 0 } @gamma),
       "positive band edge is explicit");
    my $symmetric = 1;
    for my $index (0 .. $#gamma) {
        my $opposite = $gamma[$#gamma - $index];
        $symmetric = 0
            if abs($gamma[$index][0] + $opposite->[0]) > 1e-14
            || abs($gamma[$index][1] - $opposite->[1]) > 1e-14;
    }
    ok($symmetric, "seed hybridization is particle-hole symmetric");
};

subtest "result checker accepts only a complete two-cycle transaction" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file("ITER", "2\n");
    write_file("STOP", "\n");
    write_file("nrg.calls", "1\n2\n1\n2\n");
    write_file("occupancy.log", "0 0 0 1\n0 0 0 1\n");
    write_file("param.mu.used", "0\n");
    write_file("param.mu.next", "0\n");
    write_file("DIFFS_C", "2 0.01\n");
    write_file("DIFFS_LatLoc", "2 0.02\n");
    write_file("$_", "") for qw(logNEW.1 logUNION.1 logNEW.2 logUNION.2);
    write_file(
        "OCCUPANCY_METRICS",
        "version=1\niteration=2\nmode=accurate\nstatus=within_tolerance\nerror_old=0\n" .
        "convergence_dx=0\nmu_old=0\nmu_new=0\nn_old=1\nn_new=1\n" .
        "error_new=0\ngoal=1\nweight_new=1\n",
    );
    my $gamma = "-1 0\n0 0.2\n1 0\n";
    my $real = "-1 0\n0 0\n1 0\n";
    my $imaginary = "-1 0\n0 -0.2\n1 0\n";
    for my $suffix (qw(used next)) {
        write_file("Delta.$suffix.dat", $gamma);
        write_file("ReDelta.$suffix.dat", $real);
        write_file("ImDelta.$suffix.dat", $imaginary);
        write_file("mesh.$suffix.dat", $real);
    }
    write_file("imsigma.dat", "-1 -1e-12\n0 -0.1\n1 -1e-12\n");
    write_file(
        "customfdm.avg",
        "# T Himp Hpot SdSk SigmaHd n_d\n0.01 -0.1 -0.1 -0.25 0 1\n",
    );
    for my $pair (
        ["ReDelta.dat", "ReDelta.next.dat"],
        ["ImDelta.dat", "ImDelta.next.dat"],
        ["Delta.dat", "Delta.next.dat"],
        ["param.mu", "param.mu.next"],
        ["mesh.dat", "mesh.next.dat"],
    ) {
        symlink($pair->[1], $pair->[0]) or die $!;
    }

    is(system("python3", "$fixture/check_result.py", $dir), 0,
       "complete synthetic result satisfies the checker");
    write_file("NEXT_PREPARATION_FAILED", "synthetic failure\n");
    isnt(system("python3", "$fixture/check_result.py", $dir), 0,
          "terminal publication with failed next-bath preparation is rejected");
};

chdir($original) or die $!;
done_testing();

sub read_table {
    my ($path) = @_;
    open(my $fh, "<", $path) or die "Can't read $path: $!";
    my @rows;
    while (<$fh>) {
        next if /^\s*(?:#.*)?$/;
        push(@rows, [map { 0.0 + $_ } split]);
    }
    close($fh) or die "Can't close $path: $!";
    return @rows;
}
