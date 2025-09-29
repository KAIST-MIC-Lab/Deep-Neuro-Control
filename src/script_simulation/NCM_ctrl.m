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
    mu = ncm.mu;

    m_min = 1e-5;

    alpha = ncm.alpha; 
    inv_R = ncm.inv_R;  

    if ncm.cstr_on == 0
        lbd = ncm.lbd;  
    end

    pre_W_bar = ncm.W_bar;

    sig = param.sig;
    rho = param.rho;
    beta = param.beta;

    %% INPUT MATRIX and SDC CALCULATION
    SDC = [
        -sig                     sig                0;
        rho-.5*(x(3)+xd(3))     -1                  -.5*(x(1)+xd(1));
        .5*(x(2)+xd(2))         .5*(x(1)+xd(1))     -beta
    ];
    B = param.B;

    %% OPTIMIZATION PROBLEM DEFINITION
    W_bar = sdpvar(x_num ,x_num);
    X = sdpvar(1,1);
    if ncm.cstr_on == 1
        AUX = sdpvar(x_num, x_num);

        mu = sdpvar(1,1);
        assign(mu, ncm.mu);

        assign(AUX, ncm.AUX);
    else
    % end
        mu = sdpvar(1,1);
        assign(mu, ncm.mu);

    end

    assign(W_bar, pre_W_bar);
    assign(X, ncm.X);

    % objective function
    if ncm.cstr_on == 0
        obj = X + lbd * mu;
    else
        lbd = 1e-2;
        % obj = X + trace(AUX);
        obj = X + trace(AUX)+ lbd * mu;
    end
    
    % constraints
    con = [
        -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
        eye(x_num) <= W_bar;
        W_bar <= X * eye(x_num);
        % X == mu/m_min;
        % m_min <= mu;
        1e-99 <= mu;
        1 <= X;
    ];

    if ncm.cstr_on
        max_u = param.max_u;
        con = [
            [max_u^2  (W_bar*inv(inv_R*B')*ud-mu*(x-xd))'
            (W_bar*inv(inv_R*B')*ud-mu*(x-xd)) AUX] >= 0;
            [AUX (W_bar*inv(inv_R*B'))'
            (W_bar*inv(inv_R*B')) eye(x_num)] >= 0,
            con, 
        ];
    end

    ops = sdpsettings('verbose', 0);
    ops = sdpsettings(ops, 'debug',0);

    % ops = sdpsettings(ops, 'solver','mosek');
    
    ops = sdpsettings(ops,'solver', 'sedumi');
    % ops = sdpsettings(ops,'sedumi.eps', 1e-1);
    % ops = sdpsettings(ops, 'sedumi.cg.maxiter', 1000);

    % optimize!
    sol = optimize(con, obj, ops);

    % optimization result check 
    %   (0: success, 1: infeasible, 2: unbounded, 3: max iterations, 4: numerical error)
    %   see https://yalmip.github.io/command/yalmiperror/
    if sol.problem ~= 0
    % if false
        ncm.X = ncm.X;
        ncm.W_bar = pre_W_bar;  % keep the previous value
        ncm.optDone = sol.problem;
        % warning('YALMIP Error: %s', yalmiperror);

        if ncm.cstr_on == 0
            ncm.mu = ncm.mu;
        end
    else
        ncm.X = value(X);
        ncm.W_bar = value(W_bar);
        ncm.optDone = sol.problem;

        if ncm.cstr_on == 1
            ncm.AUX = value(AUX);
            ncm.mu = value(mu);

        else
            ncm.mu = value(mu);
        end
    end

end