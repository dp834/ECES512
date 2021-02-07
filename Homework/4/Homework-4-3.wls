a = {{-.5, 0}, {0, .5}};
b = {{1},{1}};
c = {{0,1}};


contrbGram = Sum[MatrixPower[a,m].b.Transpose[b].MatrixPower[Transpose[a],m], {m,0,Infinity}]
Print["Discrete Controllability Gramian: ", contrbGram]
Print["Discrete Controllability Gramian positive definite: ", PositiveDefiniteMatrixQ[contrbGram]]

observGram = Sum[MatrixPower[Transpose[a],m].Transpose[c].c.MatrixPower[a,m], {m,0,Infinity}]
Print["Discrete Observability Gramian: ", observGram]
Print["Discrete Observability Gramian positive definite: ", PositiveDefiniteMatrixQ[observGram]]
