a = {{1,1,-2},{0,1,1},{0,0,1}}
b = {{1},{0},{1}};
c = {{2,0,0}};
k = {{k1,k2,k3}};
poles[s_] = (s-(-1 + I))(s- ( -1 - I ))(s+2) // Simplify;


siaInv = Inverse[s IdentityMatrix[3] - a] // FullSimplify
Print["Det(A):", Det[a]]
Print["(sI-A)^-1:", siaInv]
Print[TeXForm[siaInv]]

Print["A^2:", a.a]
Print[TeXForm[a.a]]

contrb = Join[b, a.b, a.a.b, 2] // Simplify;
Print["Controllability Matrix", contrb]
Print[TeXForm[contrb]]
Print["Controllability Matrix rank", MatrixRank[contrb]]

Print["Target poles:", poles[s]]

aCL = a-b.k;
aclPoly[s_] = Det[s IdentityMatrix[3] - aCL] // Simplify;

Print["(sI-(a-bk)):", Collect[aclPoly[s], s]]
Print[Simplify[poles[s] == aclPoly[s]]]
sol[s_] = Reduce[poles[s] == aclPoly[s], {k1,k2,s}];
Print["K vals:", sol[s]];

polesMatPoly = a.a.a + 4 a.a + 6 a + 4 MatrixPower[a,0]
Print["Inverse of Controllability Matrix", Inverse[contrb]]
Print["Matrix Charictaristic Poly:", polesMatPoly]
Print["C^-1 Delta(A)", Inverse[contrb].polesMatPoly]
Print["k^T:", {{0,0,1}}.Inverse[contrb].polesMatPoly]
