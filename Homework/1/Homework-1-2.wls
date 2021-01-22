gs = (3 s + 2)/(s (s^2 + 3 s + 2));
us = 1/s

Print["Partial fraction decomposition G(s) = ", Apart[gs]]
Print[TeXForm[Apart[gs]]]


Print["Partial fraction decomposition G(s)U(s) = ", Apart[gs*us]]
Print[TeXForm[Apart[gs*us]]]

yt = InverseLaplaceTransform[gs*us, s, t] //Simplify;
Print["Inverse Laplace G(s)U(s)=",  yt];
Print[TeXForm[yt]]
