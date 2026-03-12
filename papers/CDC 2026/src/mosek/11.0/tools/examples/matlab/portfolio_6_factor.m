%%
%  File : portfolio_6_factor.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 of factor type.
% 
%%
function portfolio_6_factor()

    %
    %    Purpose:
    %        Computes the optimal portfolio for a given risk 
    %     
    %    Input:
    %        n: Number of assets
    %        mu: An n dimmensional vector of expected returns
    %        G_factor_T: The factor (dense) part of the factorized risk
    %        theta: specific risk vector
    %        x0: Initial holdings 
    %        w: Initial cash holding
    %        gamma: Maximum risk (=std. dev) accepted
    %     
    %    Output:
    %        Optimal expected return and the optimal portfolio     
    %

    n      = 8;
    w      = 1.0;
    mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
    x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
    gammas = [0.24 0.28 0.32 0.36 0.4 0.44 0.48];

    B    = [ 0.4256, 0.1869; ...
             0.2413, 0.3877; ...
             0.2235, 0.3697; ...
             0.1503, 0.4612; ...
             1.5325, -0.2633; ...
             1.2741, -0.2613; ...
             0.6939, 0.2372; ...
             0.5425, 0.2116 ];

    S_F = [0.0620, 0.0577; ...
           0.0577, 0.0908 ];

    theta = [0.0720 0.0508 0.0377 0.0394 0.0663 0.0224 0.0417 0.0459];

    %% Compute the small factorization required for the model
    P  = chol(S_F)';
    G_factor  = B * P;
    G_factor_T = G_factor';

    fprintf("Basing Markowitz portfolio optimization using a factor model");
    for gamma = gammas
    
        % A model with variable x and maximized return
        model = mosekmodel(name = "Factor risk", ...
                           numvar = n, ...
                           objsense = "maximize", ...
                           objective = mu);

        % Linear budget constraint
        model.appendcons(name = "budget", F = ones(1,n), ...
                         domain = mosekdomain("equal", rhs = w + sum(x0)));

        % No shortselling, x >= 0
        model.appendcons(F = speye(n),  ...
                         domain = mosekdomain("greater than", dim = n, rhs = zeros(n,1)));

        % An affine conic constraint: [gamma, G_factor_T*x, sqrt(theta).'x ] in quadratic cone
        model.appendcons(name = "risk", ...
                         F = sparse([ zeros(1,n); G_factor_T; diag(sqrt(theta)) ]), ...
                         g = [gamma; zeros(size(G_factor_T, 1) + n, 1)], ...
                         domain = mosekdomain("qcone", dim = size(G_factor_T, 1) + n + 1 ));

        model.solve(quiet = true);

        % Check if the solution is an optimal point
        % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
        [ok,prosta,solsta] = model.hassolution("interior");
        if ~ ok || solsta ~= "OPTIMAL"
            error("Unexpected solution status");
        else
            x = model.getsolution("interior","x");
            er = mu'*x;
            fprintf("Expected return: %.4e Std. deviation: %.4e\n", er, gamma);
        end
    end
end
