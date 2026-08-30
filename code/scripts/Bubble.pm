package Bubble;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(
    bubble_options bubble_scalar chemical_potential parameter run_bubble
    write_scalar
);

my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;

sub parse_scalar {
    my ($text, $context) = @_;
    $text =~ /\A\s*($number)\s*\z/
        or die "$context did not produce exactly one numeric scalar: $text";
    my $value = 0.0 + $1;
    abs($value) <= 1.7976931348623157e308
        or die "$context produced a non-finite scalar\n";
    return $value;
}

sub parameter {
    my @names = @_;
    for my $name (@names) {
        open(my $fh, "-|", "getparam", $name, "param.loop")
            or die "Can't run getparam for $name: $!";
        local $/;
        my $output = <$fh>;
        my $success = close($fh);
        next if !$success;
        next if !defined($output) || $output =~ /\A\s*\z/;
        return parse_scalar($output, "parameter $name");
    }
    die "None of the parameters were found: " . join(", ", @names) . "\n";
}

sub chemical_potential {
    my $filename = defined($ENV{DMFT_MU_FILE}) && length($ENV{DMFT_MU_FILE})
        ? $ENV{DMFT_MU_FILE}
        : (-s "param.mu.used" ? "param.mu.used" : "param.mu");
    open(my $fh, "<", $filename) or die "Can't read $filename: $!";
    local $/;
    my $contents = <$fh>;
    close($fh) or die "Can't close $filename: $!";
    defined($contents) or die "Chemical potential not found in $filename\n";
    return parse_scalar($contents, $filename);
}

sub bubble_options {
    my %options = @_;
    my $epsabs = $options{epsabs} // "1e-9";
    my $epsrel = $options{epsrel} // "1e-8";
    my $sigma_clip = parameter("clipSigma", "clip");
    $sigma_clip > 0.0 or die "clipSigma must be positive\n";
    return (
        "--gsl-error-policy", "fail",
        "--workspace-limit", "1000",
        "--interpolation", "steffen",
        "--phi-interpolation", "steffen",
        "-k", "6",
        "-c", "30",
        "-a", $epsabs,
        "-r", $epsrel,
        "-s", $sigma_clip,
    );
}

sub bubble_scalar {
    my @arguments = @_;
    open(my $fh, "-|", "bubble", @arguments)
        or die "Can't run bubble: $!";
    local $/;
    my $output = <$fh>;
    close($fh) or die "bubble failed: $?";
    defined($output) or die "bubble produced no output\n";
    return parse_scalar($output, "bubble");
}

sub run_bubble {
    my @arguments = @_;
    system("bubble", @arguments) == 0 or die "bubble failed: $?";
}

sub write_scalar {
    my ($filename, $value) = @_;
    defined($value) && abs($value) <= 1.7976931348623157e308
        or die "Refusing to write a non-finite scalar to $filename\n";
    my $temporary = "$filename.tmp.$$";
    open(my $fh, ">", $temporary) or die "Can't write $temporary: $!";
    print {$fh} sprintf("%.17g\n", $value)
        or die "Can't write $temporary: $!";
    close($fh) or die "Can't close $temporary: $!";
    rename($temporary, $filename)
        or do {
            my $error = $!;
            unlink($temporary);
            die "Can't replace $filename: $error";
        };
}

1;
