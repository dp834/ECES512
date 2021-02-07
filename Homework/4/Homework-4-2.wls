a[t_] = {{-2, t}, {9, -3}};
b[t_] = {{0},{1}};
c[t_] = {{0,Exp[5 t]}};

Print["Controllability"]
M[0,t_] = b[t]
M[m_,t_] := -a[t].M[m-1,t] + D[M[m-1,t],t]

contrb[t_] = Join[M[0,t], M[1,t], 2]
Print["Controllability Matrix: ", contrb[t]]
Print[TeXForm[contrb[t]]]

singularT = Solve[Det[contrb[t]] == 0, {t}]
Print["Singular for values: ", singularT ]

Print["Observability"]
n[0,t_] = c[t]
n[m_,t_] := n[m-1,t].a[t] + D[n[m-1,t],t]

observ[t_] = Join[n[0,t], n[1,t], 1]
Print["Observability Matrix: ", observ[t]]
Print[TeXForm[observ[t]]]

singularT = Solve[Det[observ[t]] == 0, {t}]
Print["Singular for values: ", singularT ]
