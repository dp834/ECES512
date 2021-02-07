a = {{-1, 4}, {4, -1}};
b = {{1},{1}};
c = {{1,1}};

ctrb = Join[b, a.b, 2]
Print["Controllability Matrix: ", ctrb]
Print[TeXForm[ctrb]]

Print["Controllability Matrix Rank: ", MatrixRank[ctrb]]

pinv = {{1,1},{1,0}};
p = Inverse[pinv];
Print["P^-1: ", pinv]
Print[TeXForm[pinv]]
Print["P: ", p]
Print[TeXForm[p]]

abar = p.a.pinv;
bbar = p.b;
cbar = c.pinv;

Print["Abar: ", abar]
Print[TeXForm[abar]]
Print["Bbar: ", bbar]
Print[TeXForm[bbar]]
Print["Cbar: ", cbar]
Print[TeXForm[cbar]]
