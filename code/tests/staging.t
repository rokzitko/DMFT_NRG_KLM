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
    write_file("Delta.dat", "-1 0\n1 0\n");
    write_file("ReDelta.dat", "-1 0\n1 0\n");
    write_file("ImDelta.dat", "-1 -1\n1 -1\n");

    is(system($^X, "$scripts/next_aliases"), 0, "legacy inputs migrate");
    is(readlink("param.mu"), "param.mu.next", "param.mu is a next-input alias");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "ReDelta alias is explicit");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "ImDelta alias is explicit");
    is(readlink("Delta.dat"), "Delta.next.dat", "Gamma alias is explicit");
    ok(-s "param.mu.next", "mu target exists");
    is(system($^X, "$scripts/next_aliases"), 0, "alias maintenance is idempotent");
};

subtest "initialization treats Gamma as authoritative" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    write_file(
        "param.loop",
        "clipDelta=0.1\nbroaden_max=2\nbroaden_ratio=2\nbroaden_min=0.5\n",
    );
    write_file("param.eps", "0.25\n");
    write_file("param.mu", "0\n");
    write_file("Delta.dat", "-2 7\n-1 0.2\n1 0.3\n2 8\n");
    write_file("ReDelta.dat", "-2 9\n-1 9\n1 9\n2 9\n");
    write_file("ImDelta.dat", "-2 -9\n-1 -9\n1 -9\n2 -9\n");

    is(system($^X, "$scripts/dmft_init"), 0,
       "preexisting Gamma restart initializes successfully");
    is(readlink("Delta.dat"), "Delta.next.dat", "Gamma alias is retained");
    is(readlink("ReDelta.dat"), "ReDelta.next.dat", "real alias is regenerated");
    is(readlink("ImDelta.dat"), "ImDelta.next.dat", "imaginary alias is regenerated");
    my @gamma = read_table("Delta.next.dat");
    my @real = read_table("ReDelta.next.dat");
    my @imaginary = read_table("ImDelta.next.dat");
    is($gamma[0][1], 0, "restart projection creates the lower support guard");
    is($gamma[-1][1], 0, "restart projection creates the upper support guard");
    my $consistent = 1;
    for my $index (0 .. $#gamma) {
        $consistent = 0 if abs($gamma[$index][1] + $imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "restart ImDelta is derived from Gamma");
    ok(!grep({ abs($_->[1] - 9) < 1e-15 } @real),
       "stale restart ReDelta values are discarded");
    is(scalar(@gamma), scalar(read_table("mesh.next.dat")),
       "restart mesh follows remeshed Gamma");
};

subtest "launcher migrates an Im-only legacy result" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("res", "scripts");
    symlink("$scripts/migrate_legacy_res", "scripts/migrate_legacy_res") or die $!;
    symlink("$scripts/causalDelta", "scripts/causalDelta") or die $!;
    my $negative = "-2 -0.2\n-1 -0.2\n1 -0.2\n2 -0.2\n";
    my $positive = "-2 0.2\n-1 0.2\n1 0.2\n2 0.2\n";
    write_file("param.loop", "clipDelta=0.1\n");
    write_file("param.eps", "0.25\n");
    write_file("DOS.dat", "placeholder\n");
    write_file("CONVERGED", "\n");
    write_file("res/ImDelta.dat", $negative);
    write_file("res/$_", $positive) for qw(
        c-imG.dat c-imF.dat c-imI.dat c-reG.dat c-reF.dat c-reI.dat
        c-self.dat imsigma.dat resigma.dat imaw.dat reaw.dat
    );
    write_file("occupancy.log", "0.1 0.15 0.05 0.7\n0.15 0.2 0.05 0.8\n");

    isnt(system($^X, "$scripts/DMFT"), 0,
         "launcher reaches the convergence guard after migrating legacy results");
    ok(!-e "res", "legacy result directory is removed");
    is(read_file("param.mu.used"), "0.15\n", "used mu is recovered from occupancy history");
    my @gamma = read_table("Delta.used.dat");
    my @imaginary = read_table("ImDelta.used.dat");
    is_deeply([map { $_->[1] } @gamma], [0, 0.2, 0.2, 0],
              "legacy ImDelta is canonicalized as authoritative Gamma");
    my $consistent = 1;
    for my $index (0 .. $#gamma) {
        $consistent = 0 if abs($gamma[$index][1] + $imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "legacy ImDelta is regenerated from Gamma");
    ok(-s "ReDelta.used.dat", "causal ReDelta is generated without a legacy real part");
    is(read_file("mesh.used.dat"), read_file("Delta.used.dat"),
       "used mesh follows authoritative Gamma");
    ok(-s "c-imG.dat" && -s "self.dat", "legacy result files move to root");
};

subtest "interrupted remeshing publication" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    my $gamma = "-2 0\n-1 0.3\n1 0.3\n2 0\n";
    my $re = "-2 9\n-1 9\n1 9\n2 9\n";
    my $im = "-2 -8\n-1 -8\n1 -8\n2 -8\n";
    my $mesh = "-2 0\n-1 0\n1 0\n2 0\n";
    write_file("param.loop", "clipDelta=0.1\n");
    write_file("param.eps", "0\n");
    write_file("param.mu.next", "0\n");
    write_file("Delta.next.dat", "");
    write_file("ReDelta.next.dat", "old\n");
    write_file("ImDelta.next.dat", "old\n");
    write_file("Delta.next.dat.resampled", $gamma);
    write_file("ReDelta.next.dat.resampled", $re);
    write_file("ImDelta.next.dat.resampled", $im);
    write_file("mesh.next.dat.new", $mesh);
    write_file("mesh.next.dat.resample-ready", "ready\n");

    is(system($^X, "$scripts/dmft_init"), 0,
       "initialization recovers staged remeshing despite truncated live Gamma");
    my @restored_gamma = read_table("Delta.next.dat");
    my @restored_real = read_table("ReDelta.next.dat");
    my @restored_imaginary = read_table("ImDelta.next.dat");
    is_deeply([map { $_->[1] } @restored_gamma], [0, 0.3, 0.3, 0],
              "Gamma is restored and projected");
    ok(!grep({ abs($_->[1] - 9) < 1e-15 } @restored_real),
       "cached ReDelta is reconstructed instead of trusted");
    my $consistent = 1;
    for my $index (0 .. $#restored_gamma) {
        $consistent = 0
            if abs($restored_gamma[$index][1] + $restored_imaginary[$index][1]) > 1e-15;
    }
    ok($consistent, "cached ImDelta is reconstructed from Gamma");
    is(read_file("mesh.next.dat"), $mesh, "mesh is restored");
};

subtest "Broyden mixes Gamma with optional mu" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "res", "bin");
    my $mesh = "-2 0\n-1 0.2\n1 0.2\n2 0\n";
    my $raw = "-2 4\n-1 0.4\n1 0.4\n2 5\n";
    write_file("dmft/1-Delta.dat", $mesh);
    write_file("dmft/1-param.mu", "0\n");
    write_file("dmft/1-mu-occup.dat", "0 0.7\n");
    write_file("res/Delta.dat.OLD-common", $mesh);
    write_file("res/Delta.dat.NEW-common", $raw);
    write_file("res/param.mu.used", "0\n");
    write_file("param.eps", "0.25\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "alpha) printf 0.5;;\n" .
        "under) printf 0.5;;\n" .
        "goal) printf 0.8;;\n" .
        "maxdx) printf 0.05;;\n" .
        "broydenM) printf 50;;\n" .
        "clipDelta) printf 0.1;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--delta", "--mu",
        "--iteration", "1", "--workdir", "res",
        "--mu-input", "res/param.mu.used",
        "--mu-output", "res/param.mu.next",
        "--history-dir", "dmft", "--log", "res/occupancy.log",
    ), 0, "combined Gamma and mu update succeeds");
    my @mixed = read_table("res/Delta.dat.TEMP");
    cmp_ok(abs($mixed[1][1] - 0.3), "<=", 1e-15,
           "Broyden state contains one mixed Gamma component");
    cmp_ok(0 + read_file("res/param.mu.next"), ">", 0,
           "mu follows the Gamma block in the combined vector");

    is(system(
        $^X, "$scripts/causalDelta", "res/Delta.dat.TEMP",
        "res/Delta.next.dat", "res/ReDelta.next.dat", "res/ImDelta.next.dat",
    ), 0, "Broyden output passes through causal projection");
    my @projected = read_table("res/Delta.next.dat");
    is($projected[0][1], 0, "projection clears the lower raw endpoint");
    is($projected[-1][1], 0, "projection clears the upper raw endpoint");
    cmp_ok(abs($projected[1][1] - 0.3), "<=", 1e-15,
           "projection retains the mixed interior Gamma");
};

subtest "Broyden migrates legacy imaginary history" => sub {
    my $dir = tempdir(CLEANUP => 1);
    chdir($dir) or die $!;
    make_path("dmft", "res", "bin");
    write_file("dmft/1-Delta.dat", "-2 0\n-1 0.2\n1 0.2\n2 0\n");
    write_file("dmft/1-ImDelta.raw.dat", "-2 -0.1\n-1 -0.4\n1 -0.4\n2 -0.1\n");
    write_file("dmft/2-Delta.dat", "-2 0\n-1 0.25\n1 0.25\n2 0\n");
    write_file("res/Delta.dat.OLD-common", "-2 0\n-1 0.25\n1 0.25\n2 0\n");
    write_file("res/Delta.dat.NEW-common", "-2 0\n-1 0.5\n1 0.5\n2 0\n");
    write_file("param.loop", "");
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "alpha) printf 0.5;;\n" .
        "broydenM) printf 50;;\n" .
        "*) exit 1;;\nesac\n",
    );
    chmod(0755, "bin/getparam");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    is(system(
        "python3", "$scripts/broyden.py", "--delta", "--iteration", "2",
        "--workdir", "res", "--history-dir", "dmft",
    ), 0, "Gamma Broyden accepts converted legacy history");
    my @migrated = read_table("dmft/1-Delta.raw.dat");
    is_deeply([map { $_->[1] } @migrated], [0.1, 0.4, 0.4, 0.1],
              "legacy raw ImDelta is migrated with the correct sign");
    ok(-s "res/Delta.dat.TEMP", "migrated history contributes to a Gamma update");
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
    write_file(
        "bin/getparam",
        "#!/bin/sh\ncase \"\$1\" in\n" .
        "track) printf false;;\n" .
        "clipDelta) printf 0.1;;\n" .
        "*) exit 1;;\nesac\n",
    );
    write_file("bin/delete_NRG_output", "#!/bin/sh\nexit 0\n");
    chmod(0755, "scripts/DMFT", "bin/getparam", "bin/delete_NRG_output");
    local $ENV{PATH} = "$dir/bin:$ENV{PATH}";

    my $negative = "-2 0\n-1 -0.2\n1 -0.2\n2 0\n";
    my $positive = "-2 0\n-1 0.2\n1 0.2\n2 0\n";
    my $mesh = "-2 0\n-1 0\n1 0\n2 0\n";
    write_file("param.loop", "track=false\nclipDelta=0.1\n");
    write_file("param.eps", "0\n");
    write_causal_stage_triplet("used", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_causal_stage_triplet("next", $positive);
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

    write_file("res/ReDelta.next.dat", $mesh);
    isnt(system($^X, "$scripts/dmft_done", "--publish"), 0,
         "ready replay rejects a stale independent ReDelta");
    ok(-d "res", "rejected ready transaction remains available for repair");
    write_causal_stage_triplet("next", $positive);

    is(system($^X, "$scripts/dmft_done", "--publish"), 0, "ready stage publishes");
    ok(!-e "res", "staging directory is removed");
    is(readlink("param.mu"), "param.mu.next", "compatibility mu alias is installed");
    is(readlink("mesh.dat"), "mesh.next.dat", "compatibility mesh alias is installed");
    is(read_file("ITER"), "1\n", "iteration is published");
    is(read_file("param.mu.used"), "0.1\n", "used mu remains with results");
    is(read_file("param.mu.next"), "0.2\n", "next mu remains distinct");
    ok(-s "c-imG.dat" && -s "self.dat", "result batch is in the main directory");

    make_path("res");
    write_causal_stage_triplet("used", $positive);
    write_file("res/param.mu.used", "0.1\n");
    write_file("res/mesh.used.dat", $mesh);
    write_causal_stage_triplet("next", $positive);
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

sub write_causal_stage_triplet {
    my ($kind, $gamma) = @_;
    my $raw = "res/.Delta.$kind.raw.$$";
    write_file($raw, $gamma);
    system(
        $^X, "$scripts/causalDelta", $raw,
        "res/Delta.$kind.dat", "res/ReDelta.$kind.dat", "res/ImDelta.$kind.dat",
    ) == 0 or die "Failed to create causal $kind staging fixture: $?";
    unlink($raw) or die "Can't remove $raw: $!";
}
