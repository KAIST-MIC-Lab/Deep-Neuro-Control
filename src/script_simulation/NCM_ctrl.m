function ncm = NCM_ctrl(ncm, x, xd, ud, param)
    
    if ncm.init % for initial control (no CV-STEM)
        ncm.optDone = 0;
        ncm.init = 0;
        return

    else

    dt = ncm.dt;
    d_MAX = ncm.d_MAX;  % maximum disturbance
    m_MIN = ncm.m_MIN;

    mu = ncm.mu; 

    alpha = ncm.alpha;  % dacay rate

    inv_R = ncm.inv_R;  

    x_num = ncm.x_num;

    lbd = 2e-2; % this will be included in `ncm` structure
    bar_b1 = 3.1111e+03; % this one will be, as well

    %% 
    th = x(1);
    L = param.L;    % Inductance (mH)
    P = param.P;         % Pole pairs
    Phi = param.Phi; % Flux (Wb)
    B = [0 0; -(3*P*Phi)/(2*L)*sin(P*th) (3*P*Phi)/(2*L)*cos(P*th)]; 

    SDC = SDCmotor(x, xd, ud, param);  % System Dynamics Coefficient

    %% OPTIMIZATION PROBLEM DEFINITION
    W_bar = sdpvar(x_num ,x_num);
    X = sdpvar(1,1);
    mu = sdpvar(1,1);

    pre_W_bar = ncm.W_bar;

    assign(W_bar, pre_W_bar);
    % assign(W, eye(x_num)*0);
    assign(X, ncm.X);
    % assign(X, 1);

    obj = X + lbd * mu;
    % obj = bar_b1*d_MAX*X/alpha + lbd * mu;
    
    eps = 1e-1;
    
    con = [
        % W == W';
        -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar - eps*eye(x_num);
        % -(W_bar-eye(2))/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar - eps*eye(x_num);
        eye(x_num) <= W_bar;
        W_bar <= X * eye(x_num);
        % X == mu/m_MIN;
        1e-12 <= mu;
    ];

    ops = sdpsettings('solver', 'sedumi', ...
                  'verbose', 0, ...
                  'sedumi.eps', 1e-8);
    % optimize(con, obj);
    sol = optimize(con, obj, ops);

    if sol.problem ~= 0
        ncm.optDone = sol.problem;
        ncm.X = 0;
        ncm.mu = ncm.mu;
        ncm.W_bar = pre_W_bar;  % keep the previous value
        % warning('YALMIP Error: %s', yalmiperror);
    else
        ncm.W_bar = value(W_bar);
        ncm.X = value(X);
        ncm.mu = value(mu);
        ncm.optDone = sol.problem;
    end

    %%

end