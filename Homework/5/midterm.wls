a = {{-1, 10},{0, 1}};
b = {{-2},{0}};
c = {{-2, 3}};

Print["Eigenvalues(A): ", Eigenvalues[a]]
Print[TeXForm[Eigenvalues[a]]]

siaInv[s_] = Inverse[s * IdentityMatrix[2]-a] //Simplify
Print["(sI-A)^-1:", siaInv[s]]
Print[TeXForm[siaInv[s]]]

g[s_] = c.siaInv[s].b+d //Together

Print["g(s)", g[s]]
Print[TeXForm[g[s]]]
