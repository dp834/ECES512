a = {{1,3}, {3,1}};
b = {{1},{0}};
k = {{k1,k2}};
poles[s_] = (s-(-1 + 2 I))(s- ( -1 - 2 I )) // Simplify;

contrb = Join[b, a.b, 2] // Simplify;
Print["Controllability Matrix", contrb]
Print[TeXForm[contrb]]
Print["Controllability Matrix rank", MatrixRank[contrb]]

Print["Target poles:", poles[s]]

aCL = a-b.k;
aclPoly[s_] = Det[s IdentityMatrix[{2,2}] - aCL] // Simplify;

Print[Simplify[poles[s] == aclPoly[s]]]
sol[s_] = Reduce[poles[s] == aclPoly[s], {k1,k2,s}];
Print["K vals:", sol[s]];

polesMatPoly = a.a + 2 a + 5 MatrixPower[a,0]
Print["Matrix Charictaristic Poly:", polesMatPoly]
Print["k^T:", {{0,1}}.Inverse[contrb].polesMatPoly]
