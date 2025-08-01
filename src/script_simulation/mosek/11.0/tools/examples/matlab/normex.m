%%
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      normex.m
%
%  Purpose:   Demonstrates various norm minimization problems
%
%             * least squares regression
%             * ridge regularization
%             * lasso regularization
%             * p-norm minimization
%
%%

% This function runs a few examples
function normex()

    % Change this to display the plots
    usegraphics=false;

    % Create some random data
    N  = 100;
    x1 = sort(rand(N,1));
    p  = @(t) 1.5*t.^3-2*t.^2+0.5*t.^1-1;
    y1 = p(x1) + 0.01*randn(N,1);
    xp = smallVander(x1, 5);
    if usegraphics
      figure(1);
      scatter(x1, y1);
      hold on
    end

    % Least squares regression
    x = norm_lse(xp, y1, zeros(0,6), []);
    if usegraphics
      plot(0:0.01:1, smallVander((0:0.01:1)', 5)*x);
    end

    % With ridge regularization
    x = norm_lse_reg(xp, y1, zeros(0,6), [], 0.1);
    if usegraphics
      plot(0:0.01:1, smallVander((0:0.01:1)', 5)*x);
    end

    % Quadratic ridge
    x = norm_lse_reg_quad(xp, y1, zeros(0,6), [], 0.2);
    if usegraphics
      plot(0:0.01:1, smallVander((0:0.01:1)', 5)*x);
    end

    % Completely random large data for lasso example
    N = 100;
    K = 3000;
    F = sprandn(K, N, 0.5);
    g = randn(K, 1);
    disp(sprintf('Lasso regularization'));
    for gamma=[0.01, 0.1, 0.3, 0.6, 0.9, 1.3]
       x = norm_lse_lasso(F, g, zeros(0,N), [], gamma);
       disp(sprintf('Gamma %.4f  density %.0f%%   |Fx-g|_2: %.4f', gamma, sum(abs(x)>1e-6)/N*100, norm(F*x-g)));
    end

    % Example with the p-norm cone for various p
    % We add a far outlier to the first example
    x12 = [x1; 0.73];
    y12 = [y1; -0.99];
    xp2 = [x12.^5, x12.^4, x12.^3, x12.^2, x12.^1, x12.^0];
    parts = [];
    labels = [];

    if usegraphics
      figure(2);
      scatter(x12, y12);
      hold on
    end

    for p=[1.1, 2.0, 3.0, 6.0]
       x = norm_p_norm(xp2, y12, zeros(0,6), [], p);
       if usegraphics
         parts = [parts; plot(0:0.01:1, smallVander((0:0.01:1)', 5)*x);];
         labels = [labels; sprintf('%.1f', p)];
       end
    end
    if usegraphics
      legend(parts, labels);
    end
end

% Least squares regression
% minimize \|Fx-g\|_2
function x = norm_lse(F,g,A,b)
    n = size(F,2);
    k = size(g,1);
    m = size(A,1);

    % A model with variables x(n) and t(1)
    model = mosekmodel(name = "norm_lse", numvar = n + 1, objective = [zeros(n, 1); 1]);

    % Linear constraints
    model.appendcons(F = [A, zeros(m, 1)], domain = mosekdomain("equal", dim = m, rhs = b));

    % Affine conic constraint
    model.appendcons(F = [sparse(1,n), 1; F, sparse(k,1)], g = [0; -g], ...
                     domain = mosekdomain("quadratic cone", dim = k + 1));

    model.solve();
    xx = model.getsolution("any", "x");
    x = xx(1:n);
end

% Least squares regression with regularization
% minimize \|Fx-g\|_2 + gamma*\|x\|_2
function x = norm_lse_reg(F,g,A,b,gamma)
    n = size(F,2);
    k = size(g,1);
    m = size(A,1);

    % A model with variables x(n) and t1, t2
    model = mosekmodel(name = "norm_lse_reg", numvar = n + 2, objective = [zeros(n, 1); 1; gamma]);

    % Linear constraints
    model.appendcons(F = [A, zeros(m, 2)], domain = mosekdomain("equal", dim = m, rhs = b));

    % Affine conic constraint for \|Fx-g\|_2
    model.appendcons(F = [sparse(1,n),        1, 0; ...
                          F,            sparse(k,2)], ...
                     g = [0; -g], ...
                     domain = mosekdomain("quadratic cone", dim = k + 1));

    % Affine conic constraint for \|x\|_2
    model.appendcons(F = [sparse(1, n+1),        1; ...
                          speye(n),     sparse(n,2) ], ...
                     domain = mosekdomain("quadratic cone", dim = n + 1));

    model.solve();
    xx = model.getsolution("any", "x");
    x = xx(1:n);
end

% Least squares regression with regularization
% The "classical" quadratic version
% minimize \|Fx-g\|_2^2 + gamma*\|x\|_2^2
function x = norm_lse_reg_quad(F,g,A,b,gamma)
    n = size(F,2);
    k = size(g,1);
    m = size(A,1);

    % A model with variables x(n) and t1, t2
    model = mosekmodel(name = "norm_lse_reg_quad", numvar = n + 2, objective = [zeros(n, 1); 1; gamma]);

    % Linear constraints
    model.appendcons(F = [A, zeros(m, 2)], domain = mosekdomain("equal", dim = m, rhs = b));

    % Affine conic constraint for \|Fx-g\|_2^2
    model.appendcons(F = [sparse(1,n) 1 0; ...
                          sparse(1, n+2); ...
                          F sparse(k, 2)], ...
                     g = [0; 0.5; -g], ...
                     domain = mosekdomain("rotated quadratic cone", dim = k + 2));

    % Affine conic constraint for \|x\|_2^2
    model.appendcons(F = [sparse(1, n+1) 1; ...
                          sparse(1, n+2); ...
                          speye(n) sparse(n, 2) ], ...
                     g = [0; 0.5; zeros(n, 1)], ...
                     domain = mosekdomain("rotated quadratic cone", dim = n + 2));

    model.solve();
    xx = model.getsolution("any", "x");
    x = xx(1:n);
end

% Least squares regression with lasso regularization
% minimize \|Fx-g\|_2 + gamma*\|x\|_1
function x = norm_lse_lasso(F,g,A,b,gamma)
    n = size(F,2);
    k = size(g,1);
    m = size(A,1);

    % A model with variables x(n), u(n) and t1, t2
    model = mosekmodel(name = "norm_lse_lasso", numvar = n + n + 2, objective = [zeros(2*n, 1); 1; gamma]);

    % Linear constraints
    model.appendcons(F = [A, zeros(m, 2)], domain = mosekdomain("equal", dim = m, rhs = b));

    % u >=. |x|
    model.appendcons(F = [speye(n)  speye(n) sparse(n, 2); ...
                          -speye(n) speye(n) sparse(n, 2)], ...
                     domain = mosekdomain("greater than", dim = 2*n, rhs = zeros(2*n, 1)));

    % t2 >= sum(u)
    model.appendcons(F = [sparse(1, n), -ones(1, n), 0, 1], ...
                     domain = mosekdomain("greater than", rhs = 0));

    % Affine conic constraint for \|Fx-g\|_2
    model.appendcons(F = [sparse(1, 2*n) 1 0; ...
                          F sparse(k, n + 2)], ...
                     g = [0; -g], ...
                     domain = mosekdomain("quadratic cone", dim = k + 1));

    model.solve();
    xx = model.getsolution("any", "x");
    x = xx(1:n);
end

% P-norm minimization
% minimize \|Fx-g\|_p
function x = norm_p_norm(F,g,A,b,p)
    n = size(F,2);
    k = size(g,1);
    m = size(A,1);

    % A model with variables x(n), r(k) and t(1)
    model = mosekmodel(name = "norm_p_norm", numvar = n + k + 1, objective = [zeros(n + k, 1); 1]);

    % Linear constraints
    model.appendcons(F = [A, sparse(m, k + 1)], domain = mosekdomain("equal", dim = m, rhs = b));

    % t == sum(r)
    model.appendcons(F = [sparse(1, n), ones(1, k), -1], ...
                     domain = mosekdomain("equal", rhs = 0));

    % Permutation matrix which picks triples (r_i, t, F_ix-g_i)
    M = [];
    for i=1:3
      M = [M, sparse(i:3:3*k, 1:k, ones(k,1), 3*k, k)];
    end

    Fcon = M * [sparse(k, n), speye(k), sparse(k,1); ...
                sparse(k, n+k), ones(k, 1); ...
                F, sparse(k, k+1)];
    gcon = M * [sparse(2*k, 1); -g];

    model.appendcons(F = Fcon, g = gcon, ...
                     domain = mosekdomain("power cone", n = k, dim = 3, alpha = [1; p-1]));

    model.solve();
    xx = model.getsolution("any", "x");
    x = xx(1:n);
end

% Just a few columns of the Vandermode matrix
function M = smallVander(v, k)
M = [];
for i=0:k 
    M = [v.^i, M];
end
end
