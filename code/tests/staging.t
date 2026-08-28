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
    write_file("ReDelta.dat", "-1 0\n1 0\n");
    write_file("ImDelta.dat", "-1 -1\n1 -1\n");

    is(system($^X, "$scripts/next_aliases"), 0, "legacy inputs migrate");
    is(readlink("param.mu"), "param.mu.next", "param.mu is a next-input alias");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "ReDelta alias is explicit");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "ImDelta alias is explicit");
    ok(-s "param.mu.next", "mu target exists");
    is(system($^X, "$scripts/next_aliases"), 0, "alias maintenance is idempotent");
};

subtest "legacy result migration" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("res");
    my $negative = "-1 -0.2\n1 -0.2\n";
    my $positive = "-1 0.2\n1 0.2\n";
    write_file("res/ReDelta.dat", "-1 0\n1 0\n");
    write_file("res/ImDelta.dat", $negative);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("occupancy.log", "0.1 0.15 0.05 0.7\n0.15 0.2 0.05 0.8\n");

    is(system($^X, "$scripts/migrate_legacy_res"), 0, "legacy res batch migrates");
    ok(!-e "res", "legacy result directory is removed");
    is(read_file("param.mu.used"), "0.15\n", "used mu is recovered from occupancy history");
    is(read_file("Delta.used.dat"), $positive, "used solver Delta is derived from ImDelta");
    ok(-s "c-imG.dat" && -s "self.dat", "legacy result files move to root");
};

subtest "interrupted remeshing publication" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    my $re = "-1 0.1\n1 0.1\n";
    my $im = "-1 -0.3\n1 -0.3\n";
    my $mesh = "-1 0\n1 0\n";
    write_file("ReDelta.next.dat", "");
    write_file("ImDelta.next.dat", "-1 -9\n1 -9\n");
    write_file("ReDelta.next.dat.resampled", $re);
    write_file("ImDelta.next.dat.resampled", $im);
    write_file("mesh.next.dat.new", $mesh);
    write_file("mesh.next.dat.resample-ready", "ready\n");

    is(system(
        $^X, "$scripts/resampleDelta",
        "ReDelta.next.dat", "ImDelta.next.dat", "mesh.next.dat",
    ), 0, "completed remeshing stage repairs a truncated live input");
    is(read_file("ReDelta.next.dat"), $re, "ReDelta is restored");
    is(read_file("ImDelta.next.dat"), $im, "ImDelta is restored");
    is(read_file("mesh.next.dat"), $mesh, "mesh is restored");
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

subtest "ready publication and cleanup" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("res", "scripts", "bin");
    symlink("$scripts/copyresults", "scripts/copyresults") or die $!;
    symlink("$scripts/next_aliases", "scripts/next_aliases") or die $!;

    write_file("scripts/DMFT", "#!/bin/sh\nexit 0\n");
    write_file("bin/getparam", "#!/bin/sh\nprintf false\n");
    write_file("bin/delete_NRG_output", "#!/bin/sh\nexit 0\n");
    chmod(0755, "scripts/DMFT", "bin/getparam", "bin/delete_NRG_output");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    my $negative = "-1 -0.2\n1 -0.2\n";
    my $positive = "-1 0.2\n1 0.2\n";
    my $mesh = "-1 0\n1 0\n";
    write_file("res/ReDelta.used.dat", $mesh);
    write_file("res/ImDelta.used.dat", $negative);
    write_file("res/Delta.used.dat", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_file("res/ReDelta.next.dat", $mesh);
    write_file("res/ImDelta.next.dat", $negative);
    write_file("res/Delta.next.dat", $positive);
    write_file("res/param.mu.next", "0.2\n");
    write_file("res/mesh.next.dat", $mesh);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("res/$_", "value\n") for qw(custom custom.avg customfdm customfdm.avg);
    write_file("res/ITER", "1\n");
    write_file("res/DECISION", "continue\n");
    write_file("res/READY", "source_iter=0\ntarget_iter=1\ndecision=continue\n");
    write_file("param.loop", "track=false\n");

    is(system($^X, "$scripts/dmft_done", "--publish"), 0, "ready stage publishes");
    ok(!-e "res", "staging directory is removed");
    is(readlink("param.mu"), "param.mu.next", "compatibility mu alias is installed");
    is(readlink("mesh.dat"), "mesh.next.dat", "compatibility mesh alias is installed");
    is(read_file("ITER"), "1\n", "iteration is published");
    is(read_file("param.mu.used"), "0.1\n", "used mu remains with results");
    is(read_file("param.mu.next"), "0.2\n", "next mu remains distinct");
    ok(-s "c-imG.dat" && -s "self.dat", "result batch is in the main directory");

    make_path("res");
    write_file("res/ReDelta.used.dat", $mesh);
    write_file("res/ImDelta.used.dat", $negative);
    write_file("res/Delta.used.dat", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_file("res/ReDelta.next.dat", $mesh);
    write_file("res/ImDelta.next.dat", $negative);
    write_file("res/Delta.next.dat", $positive);
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
