%%
%  File : affco2.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :
%    Implements a basic tutorial example with affine conic constraints:
%
%    minimize t
%    st.      (d, z1*y1,... zn*yn) \in Q^{n+1}
%             (yi, 1, ai*t)        \in EXP, i=1,\ldots,n
%
%    with input ai<0, zi, d.
%
%    See also https://docs.mosek.com/modeling-cookbook/expo.html#hitting-time-of-a-linear-system
% 
%%
function affco2()

n = 3;
z = [3.3, 2.2, 1.3]';
a = [-0.08, -0.3, -0.06]';
d = 0.5;

t = firstHittingTime(n, z, a, d);
disp(sprintf('\nt = %.4e', t))


function t = firstHittingTime(n, z, a, d)

   % Variables [t, y1, ..., yn]
   model = mosekmodel(...
       name = "affco2", ...
       numvar = n + 1, ...
       objective = [1 zeros(1,n)], ...
       objsense = "minimize");

   model.varname([1], ["t"]);

   % Quadratic cone
   model.appendcons(F = diag([0; z]), ...
                    g = [d; zeros(n,1)], ...
                    domain = mosekdomain("qcone", dim = n + 1));

   % All exponential cones (their number is n)
   FExp = sparse([1:3:3*n    3:3:3*n], ...
                 [2:n+1      ones(1,n)], ...
                 [ones(1,n)  a']);
   gExp = repmat([0; 1; 0], n, 1);
   model.appendcons(F = FExp, g = gExp, ...
                    domain = mosekdomain("exp", dim = 3, n = n));

   % Solve and get solution
   model.solve();

   if model.hassolution("interior")
      [x, prosta, solsta] = model.getsolution("interior", "x");
      t = x(1);
   end
