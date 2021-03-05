a = {{-100, 0,0},{0,0,1},{-2.803,982,0}};
b = {{1/.01}, {0}, {0}}
c = {{0,1,0}};
d = {{0}};
x0 = {{7},{.0098},{0}}

contrPoles[s_] = (s-(-20+20 I))(s-(-20-20 I))(s-(-100)) //Expand

Print["Eigenvalues a:", Eigenvalues[a]]

csiab = c.Inverse[s IdentityMatrix[{3,3}] - a].b // Factor
Print["Transfer Function:", csiab]
Print[TeXForm[csiab]]

contrb = Join[b, Join[a.b, a.a.b,2], 2] // Simplify;
Print["Controllability Matrix", contrb]
Print[TeXForm[contrb]]
Print["Controllability Matrix rank", MatrixRank[contrb]]

observ = Join[c, Join[c.a, c.a.a, 1] , 1] // Simplify;
Print["Observability Matrix:,", observ];
Print[TeXForm[observ]];
Print["Rank(Observability Matrix): ", MatrixRank[observ]];

Print["Controlability poles: ", contrPoles[s]]

contrPolePolyMatrix = a.a.a + 140 a.a + 4800 a + 80000 IdentityMatrix[3];
kvals = {{0,0,1}}.Inverse[contrb].contrPolePolyMatrix;
Print["Gains matrix k: ", kvals]
Print[TeXForm[kvals]]

closedLoopResponse = c.Inverse[s IdentityMatrix[3] -(a-b.kvals)].b //Simplify;
Print["Closed Loop Response: ", closedLoopResponse];
Print[TeXForm[closedLoopResponse]];

obsrvPoles[s_] = (s-5(-20+20 I))(s-5(-20-20 I))(s-5(-100)) // Expand

Print["Observability poles: ", obsrvPoles[s]]
observPolePolyMatrix = a.a.a + 700 a.a + 120000 a + 10000000 IdentityMatrix[3];
lvals = observPolePolyMatrix.Inverse[observ].{{0},{0},{1}};
Print["Gains matrix l: ", lvals]
Print[TeXForm[lvals]]
