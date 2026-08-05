def1ch[1];

snegrealconstants[eps1, U1, B, JK];

If[!paramexists["spin", "extra"],
   MyError["Define the spin of the impurity!"];
];
SPIN = ToExpression @ param["spin", "extra"];
MyPrint["SPIN=", SPIN];

Module[{sz, sp, sm, sx, sy, oz, op, om, ss},
       sz = spinketbraZ[SPIN];
       sp = spinketbraP[SPIN];
       sm = spinketbraM[SPIN];
       sx = spinketbraX[SPIN];
       sy = spinketbraY[SPIN];

       oz = nc[ sz, spinz[ d[] ] ];
       op = nc[ sp, spinminus[ d[] ] ];
       om = nc[ sm, spinplus[ d[] ] ];

       ss = oz + 1/2 (op + om) // Expand;

       Heps = eps1 number[d[]] + B spinz[d[]] + B sz;
       Hpot = U1 hubbard[d[]] + JK ss; (* No field! *)
       Himp = Heps + Hpot; (* Includes field and eps *)

       H = H0 + Himp + Hc;
];

MAKESPINKET = SPIN;

(* All operators which contain d[], except hybridization (Hc). *)
Hselfd = Hpot;

selfopd = ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];

SigmaHd = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #1] ] /. params &;
SigmaHdAvg := Expand @ (SigmaHd[UP]+SigmaHd[DO])/2;

Print["SigmaH[UP]=", SigmaHd[UP] ];
Print["SigmaH[DO]=", SigmaHd[DO] ];
Print["SigmaH=", SigmaHdAvg ];
