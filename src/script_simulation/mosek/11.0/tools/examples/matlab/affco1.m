%%
%  File : affco1.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :
%    Implements a basic tutorial example with affine conic constraints:
%
%    maximize x_1^(1/3) + (x_1+x_2+0.1)^(1/4)
%    st.      (x_1-0.5)^2 + (x_2-0.6)^2 <= 1
%             0 <= x_1 <= x_2 + 1
% 
%%
function affco1()

   % Variables [x1; x2; t1; t2]
   model = mosekmodel(name = "affco1", numvar = 4, ...
                      objsense = "maximize", objective = [0, 0, 1, 1], ...
                      F = sparse([zeros(1,4); speye(2) zeros(2,2)]), ...    % The quadratic cone constraint
                      g = [1 -0.5 -0.6]', ...
                      domain = mosekdomain("qcone", dim = 3));

   % The power cones added as one block:
   model.appendcons(name="pow", ...
                    F = sparse([1,3,4,4,6], [1,3,1,2,4], ones(1,5)), ...
                    g = [0 1 0 0.1 1 0]', ...
                    domain = [mosekdomain("pow", dim = 3, alpha = [1 2]'), ...     % Exponents [ 1/3, 2/3 ]
                              mosekdomain("pow", dim = 3, alpha = 0.25) ]);        % Exponents [ 0.25, 0.75 ]

   % Linear inequality x_1 - x_2 <= 1
   model.appendcons(F = [1 -1 0 0], domain = mosekdomain("less than", rhs = 1));

   % Solve and get solution
   model.solve();

   if model.hassolution("interior")
      [xx, prosta, solsta] = model.getsolution("interior", "x");
      x = xx(1:2);
      t = xx(3:4);
      disp(x);
   end
end
