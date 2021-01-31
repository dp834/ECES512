eqns = {x1'[t] == -x1[t] + Exp[2 t] x2[t], x2'[t] == -x2[t]}
Print[eqns]
difeqSol = DSolve[eqns, {x1, x2}, t] // Simplify
Print["Differential Equation Solutions:", difeqSol]
Print[TeXForm[difeqSol]]

a[t_] = {{Exp[-t], Exp[t]/2},{0, Exp[-t]}};
Print["X(t):", a[t]]
Print["X^-1(t):",Simplify[Inverse[a[t]]]]
Print["Phi(t,t0):",Simplify[a[t].Inverse[a[t0]]]]
Print[TeXForm[Simplify[a[t].Inverse[a[t0]]]]]
