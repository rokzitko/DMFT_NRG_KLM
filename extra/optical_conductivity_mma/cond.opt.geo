#!/usr/bin/env perl
# Rok Zitko, 2016-2026
# Bethe lattice
# This version supports shift by param.eps and param.mu.
# Results in units of sigma_MIR.

chomp(my $TD = `getparam T param.loop`);
print "T/D=$TD\n";
$beta = 1/$TD;

# Chemical potential
my $mu_file = -s "param.mu.used" ? "param.mu.used" : "param.mu";
my $mu = `cat $mu_file`;
chomp($mu);

# epsilon_d
my $eps = `cat param.eps`;
chomp($eps);

my $shift = $mu - $eps;
print "mu=$mu eps=$eps shift=$shift\n";

open(F, ">tmp.tmp") or die;
print F <<ENDSCRIPT;

Off[NIntegrate::slwcon];

Phi[eps_] = (1-eps^2)^(3/2);

CLIP = 10^-16;
clip[x_] := If[x > -CLIP, -CLIP, x];

t = Import["imsigma.dat"];
t = Map[{First[#], clip[Last[#]]}&, t];
im = Map[{First[#], I Last[#]}&, t];
t = Import["resigma.dat"];
re = Map[{First[#],   Last[#]}&, t];
t = re + im;
t[[All,1]] = re[[All,1]];

\\[CapitalSigma] = Interpolation[t, InterpolationOrder -> 2];

\\[Beta] = $beta;
T=1/\\[Beta];
shift = $shift;

f[x_] := 1/(1+Exp[\\[Beta] x]);

calc[o_] := Module[{},
 F = (f[\\[Omega]]-f[\\[Omega]+o])/o;

 cond = NIntegrate[ F Phi[\\[Epsilon]] *
    Im[1/( \\[Omega] + shift   - \\[Epsilon] - \\[CapitalSigma] [\\[Omega]] )] *
    Im[1/((\\[Omega]+o) + shift - \\[Epsilon] - \\[CapitalSigma] [\\[Omega]+o] )],
          {\\[Omega], -20T-o, -o, 0, 20T}, {\\[Epsilon], -1, 0, 1}, MaxRecursion -> 12,
  WorkingPrecision -> 14];

 condMIR = cond * 2/Pi
];

r=1.5;
imax = Floor[20/Log[r]];

omlist = Table[ 10.^-8 r^i, {i, 0, imax}];

tab = Map[Echo[{#, calc[#]}] &, omlist];
Export["cond.opt.geo.dat", tab];

ENDSCRIPT

close(F);

system "math -batchinput -batchoutput <tmp.tmp";
unlink "tmp.tmp";
