#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Find qw(find);
use File::Spec;
use FindBin;
use Test::More;

my $root = abs_path("$FindBin::Bin/../..");
my %files = (
    perl => [],
    python => [],
    shell => [],
);

find(
    {
        no_chdir => 1,
        wanted => sub {
            my $path = $File::Find::name;
            if (-d $path) {
                $File::Find::prune = 1
                    if $path =~ m{/(?:\.git|__pycache__)\z};
                return;
            }
            return if !-f $path;

            open(my $fh, "<", $path) or die "Can't read $path: $!";
            my $first = <$fh> // "";
            close($fh) or die "Can't close $path: $!";

            if ($path =~ /\.(?:pm|t)\z/ || $first =~ /^\#!.*\bperl(?:\s|$)/) {
                push(@{$files{perl}}, $path);
            } elsif ($path =~ /\.py\z/ || $first =~ /^\#!.*\bpython3?(?:\s|$)/) {
                push(@{$files{python}}, $path);
            } elsif ($first =~ m{^\#!.*(?:/|\s)(?:ba)?sh(?:\s|$)}) {
                push(@{$files{shell}}, $path);
            }
        },
    },
    map { "$root/$_" } qw(code extra more_examples),
);

subtest "Perl syntax" => sub {
    ok(@{$files{perl}} > 0, "found Perl sources");
    for my $path (sort @{$files{perl}}) {
        is(system($^X, "-c", $path), 0, relative_path($path));
    }
};

subtest "Python syntax" => sub {
    ok(@{$files{python}} > 0, "found Python sources");
    my $compile = <<'PYTHON';
import pathlib
import sys
compile(pathlib.Path(sys.argv[1]).read_bytes(), sys.argv[1], "exec")
PYTHON
    for my $path (sort @{$files{python}}) {
        is(system("python3", "-c", $compile, $path), 0, relative_path($path));
    }
};

subtest "shell syntax" => sub {
    ok(@{$files{shell}} > 0, "found shell sources");
    for my $path (sort @{$files{shell}}) {
        is(system("bash", "-n", $path), 0, relative_path($path));
    }
};

done_testing();

sub relative_path {
    return File::Spec->abs2rel($_[0], $root);
}
