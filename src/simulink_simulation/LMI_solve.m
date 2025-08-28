function [W,X,optDone] = LMI_solve(pre_W, pre_X, SDC, B, inv_R, dt)
    
    %% 
    x_num = 3;

    % m_MAX = 0e3;  % maximum mass
    % m_MIN = 0e-1;   % minimum mass
    d_MAX = 1;  % maximum disturbance
    
    mu = 1e0;  
    alpha = 1e-1;  % decay rate
    
    %% OPTIMIZATION PROBLEM DEFINITION
    W = sdpvar(x_num ,x_num);
    assign(W, pre_W);
    X = sdpvar(1,1);
    % mu = sdpvar(1,1);
    assign(X, pre_X);
    % assign(mu, 1);

    c0 = 1;
    c1 = 0;

    obj = c0*X;
    % obj = c0*X + c1*mu;
    
    % if ncm.init
    if false
        con = [
                % W == W';
                (W'*SDC' + SDC*W) - 2*d_MAX*B*inv_R*B' <= -2*alpha*W;
                eye(x_num) <= mu*W;
                mu*W <= X * eye(x_num);
                % X == mu/m_MIN;
            ];
    else
        con = [
            % W == W';
            -(W-pre_W)/dt + (W'*SDC' + SDC*W) - 2*d_MAX*B*inv_R*B' <= -2*alpha*W;
            eye(x_num) <= mu*W;
            mu*W <= X * eye(x_num);
            % X == mu/m_MIN;
        ];

        init = 0;
    end

    ops = sdpsettings('solver', 'sedumi', ...
                  'verbose', 0, ...
                  'sedumi.eps', 1e-8);
    % optimize(con, obj);
    sol = optimize(con, obj, ops);

    if sol.problem ~= 0
        optDone = 0;
        X = 0;
        W = pre_W;  % keep the previous value
        % warning('YALMIP Error: %s', yalmiperror);
    else
        W = value(W);
        X = value(X);
        optDone = 1;
    end

end
