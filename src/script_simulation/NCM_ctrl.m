function ncm = NCM_ctrl(ncm, x, xd, ud, param)
    
    % for initial control (no CV-STEM)
    if ncm.init 
    % if true
        ncm.optDone = 0;
        ncm.init = 0;
        return
    end

    %% FROM STRUCTURE TO LOCAL VARIABLES
    dt = ncm.dt;
    x_num = ncm.x_num;
    
    alpha = ncm.alpha; 
    d_MAX = ncm.d_MAX; 

    inv_R = ncm.inv_R;  

    lbd = ncm.lbd;  

    %% INPUT MATRIX and SDC CALCULATION
    th = x(1);
    J = param.J;    
    P = param.P;     
    Phi = param.Phi;
    B = [0 0; -(3*P*Phi)/(2*J)*sin(P*th) (3*P*Phi)/(2*J)*cos(P*th)]; 

    SDC = SDCmotor(x, xd, ud, param);  

    %% OPTIMIZATION PROBLEM DEFINITION
    W_bar = sdpvar(x_num ,x_num);
    X = sdpvar(1,1);
    mu = sdpvar(1,1);

    pre_W_bar = ncm.W_bar;

    assign(W_bar, pre_W_bar);
    assign(X, ncm.X);
    assign(mu, ncm.mu);
    
    % objective function
    obj = X + lbd * mu;
    
    % numerical trouble-shooting
    eps = 0e-1;
    m_min = 1e-12;

    % constraints
    con = [
        -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar - eps*eye(x_num);
        % -(W_bar-eye(x_num))/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar - eps*eye(x_num);
        eye(x_num) <= W_bar;
        W_bar <= X * eye(x_num);
        % X == mu/m_min;
        m_min <= mu;
    ];

    ops = sdpsettings('verbose', 0);
    ops = sdpsettings(ops, 'debug',0);

    % ops = sdpsettings(ops, 'solver','mosek');
    
    ops = sdpsettings(ops,'solver', 'sedumi');
    ops = sdpsettings(ops,'sedumi.eps', 1e-8);
    ops = sdpsettings(ops, 'sedumi.cg.maxiter', 1000);

    % optimize!
    sol = optimize(con, obj, ops);

    % optimization result check 
    %   (0: success, 1: infeasible, 2: unbounded, 3: max iterations, 4: numerical error)
    %   see https://yalmip.github.io/command/yalmiperror/
    if sol.problem ~= 0
    % if false
        ncm.X = ncm.X;
        ncm.mu = ncm.mu;
        ncm.W_bar = pre_W_bar;  % keep the previous value
        ncm.optDone = sol.problem;
        % warning('YALMIP Error: %s', yalmiperror);
    else
        ncm.X = value(X);
        ncm.mu = value(mu);
        ncm.W_bar = value(W_bar);
        ncm.optDone = sol.problem;
    end

end