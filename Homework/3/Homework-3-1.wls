a = {{0, 1, 2},{0, -2, 3},{0, 0, 0}};
b = {{1},{1},{0}};
c = {{0, 1, 0}};
d = {{2}};

Print["Eigenvalues(A): ", Eigenvalues[a]]
Print[TeXForm[Eigenvalues[a]]]

siaInv[s_] = Inverse[s * IdentityMatrix[3]-a] //Simplify
Print["(sI-A)^-1:", siaInv[s]]
Print[TeXForm[siaInv[s]]]

g[s_] = c.siaInv[s].b+d //Together

Print["g(s)", g[s]]
Print[TeXForm[g[s]]]
