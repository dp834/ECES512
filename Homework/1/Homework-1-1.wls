a = {{-2,0},{1,-1}};
b = {{1},{0}};

eat = MatrixExp[a t] //FullSimplify;
Print["e^(At): ", eat];
Print[TeXForm[eat]];
toIntegrate = MatrixExp[a (t-tau)].b // FullSimplify;
Print["e^(A(t-tau))B: ", toIntegrate];
Print[TeXForm[toIntegrate]]

int = Integrate[toIntegrate, {tau, 0, t}]
Print["Convolution e^(A(t-tau)) B: ", int];
Print[TeXForm[int]]

xT = eat.{{2},{3}} + int //Expand;
Print["x(t) = ", xT];
Print[TeXForm[xT]]
