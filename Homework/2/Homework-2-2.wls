A = {{-100, 0, 0}, {0, 0, 1}, {-2.803, 982, 0}};
B = {{1/.01}, {0}, {0}};
c = {0, 1, 0};


MatrixExp[A (t - tau)].B 7 // FullSimplify // MatrixForm
intsol = Integrate[MatrixExp[A (t - tau)].B 7, {tau, 0, t}] // Chop // FullSimplify
Print["Convolution Solution:", intsol]
Print[TeXForm[intsol]]

x0 = {{7},{.01},{0}};
Assumptions -> t >= 0;
xt = MatrixExp[A t].x0 + intsol // Expand // Chop // FullSimplify
Print["x(t):", xt]
Print[TeXForm[xt]]

Print["Matlab Form:"]
<< "ToMatlab`"
Print[ToMatlab[xt]]
