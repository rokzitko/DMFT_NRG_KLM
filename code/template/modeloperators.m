Module[{t},
  t = {};
  MPVCFAST = False;
  t = Join[t, mtSingletOp["SigmaHd",   SigmaHdAvg ]];
  t = Join[t, mtSingletOp["SigmaHd-u", SigmaHd[UP] ]];
  t = Join[t, mtSingletOp["SigmaHd-d", SigmaHd[DO] ]];
  MPVCFAST = True;
  (* Spin correlation between d orbital and spin *)
  If[calcopq["SdSk"],
    Module[{SPIN, sx, sy, sz, op},
      SPIN = ToExpression @ param["spin", "extra"];
      sx = spinketbraX[SPIN];
      sy = spinketbraY[SPIN];
      sz = spinketbraZ[SPIN];
      op = sx~nc~spinx[d[]] + sy~nc~spiny[d[]] + sz~nc~spinz[d[]];
      t = Join[t, mtSingletOp["SdSk", op] ];
    ];
  ];
  texportable = t;
];
texportable
