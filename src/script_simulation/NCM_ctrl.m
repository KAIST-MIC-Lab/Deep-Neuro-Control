    function ncm = NCM_ctrl(ncm, x, xd, ud, param)
    % for initial control
    if ncm.init 
        ncm.optDone = 0;
        ncm.init = 0;
        ncm.contraction_flag = 1;
        return
    end

    %% FROM STRUCTURE TO LOCAL VARIABLES
    dt = ncm.dt;
    x_num = ncm.x_num;

    m_min = 1e-5;

    alpha = ncm.alpha; 
    inv_R = ncm.inv_R;  

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
    switch ncm.ctrl_no 
        case 0 % existing
            pre_W_bar = ncm.W_bar; pre_mu = ncm.mu;

            W_bar = sdpvar(x_num ,x_num); assign(W_bar, pre_W_bar);
            mu = sdpvar(1,1); assign(mu, pre_mu);
            X = sdpvar(1,1); assign(X, ncm.X);

            % objective function
            obj = X + ncm.lbd * mu;
        
            % constraint
            con = [
                -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
                eye(x_num) <= W_bar;
                W_bar <= X * eye(x_num);
                % X == mu/m_min;
                % m_min <= mu;
                1e-99 <= mu;
                1 <= X;
            ];

        case 1 % proposed 1 (effective space)
            pre_M = ncm.M; pre_M_bar = ncm.M_bar;
            e = x-xd;

            M = sdpvar(x_num, x_num); assign(M, pre_M);
            M_bar = sdpvar(x_num, x_num); assign(M_bar, pre_M_bar);
            mu = sdpvar(1,1); assign(mu, ncm.mu);

            y = sdpvar(x_num,1); assign(y, pre_M*e);

            mu_min = sdpvar(1,1);
            mu_max = sdpvar(1,1);

            % objective function
            obj = 1;
            obj = -mu;
            obj = norm(M-eye(x_num));
            % obj = norm(M-eye(x_num)) - mu;
            % obj = norm(M-eye(x_num)*mu) - mu;
            % obj = norm(M_bar-eye(x_num)) - mu;
            % obj = norm(M - (trace(M)/x_num)*eye(x_num)) - (trace(M)/x_num) * 1e5; 
            % obj = norm(M - (trace(M)/x_num)*eye(x_num), 'fro')^2 - (trace(M)/x_num) * 1e1;

            % obj = trace(M) - log(det(M));
            % obj = log(mu_max) - log(mu_min);

            % constraint
            u_max = param.max_u;
            con = [
                -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num))*y - (e'*pre_M*e)/dt <= 0;
                y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
                M*e == y;
                % M >= mu * eye(x_num);
                % mu >= 1e-99
                % M <= 1e9 * eye(x_num);
                % M == mu * M_bar;
            ];

        case 2 % proposed 2 (inverse)
            % AUX = sdpvar(x_num, x_num);
    
            % mu = sdpvar(1,1);
            % assign(mu, ncm.mu);
    
            % assign(AUX, ncm.AUX);

            % % objective function
            % lbd = 1e-2;
            % % obj = X + trace(AUX);
            % obj = X + trace(AUX)+ lbd * mu;

            % % constraints
            % con = [
            %     -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
            %     eye(x_num) <= W_bar;
            %     W_bar <= X * eye(x_num);
            %     % X == mu/m_min;
            %     % m_min <= mu;
            %     1e-99 <= mu;
            %     1 <= X;
            % ];
        
            % max_u = param.max_u;
            % con = [
            %     [max_u^2  (W_bar*inv(inv_R*B')*ud-mu*(x-xd))'
            %     (W_bar*inv(inv_R*B')*ud-mu*(x-xd)) AUX] >= 0;
            %     [AUX (W_bar*inv(inv_R*B'))'
            %     (W_bar*inv(inv_R*B')) eye(x_num)] >= 0,
            %     con, 
            % ];
    end
    
    %% OPTIMIZATION 
    switch ncm.ctrl_no
        case 0
            ops = sdpsettings('verbose', 0);
            ops = sdpsettings(ops, 'debug',0);
    
            % ops = sdpsettings(ops, 'solver','mosek');
            ops = sdpsettings(ops,'solver', 'sedumi');
            % ops = sdpsettings(ops,'sedumi.eps', 1e-1);
            % ops = sdpsettings(ops, 'sedumi.cg.maxiter', 1000);
            
        case 1 
            ops = sdpsettings('solver', 'fmincon', 'verbose', 0, 'fmincon.Algorithm', 'sqp');
            % ops = sdpsettings('verbose', 0);
            % ops = sdpsettings('solver', 'ipopt');

    end

    % optimize!
    sol = optimize(con, obj, ops);

    % optimization result check 
    %   (0: success, 1: infeasible, 2: unbounded, 3: max iterations, 4: numerical error)
    %   see https://yalmip.github.io/command/yalmiperror/

    %% TERMINIATING OPTIMIZATION
    if sol.problem ~= 0
        ncm.optDone = sol.problem;
        return;
    end

    switch ncm.ctrl_no
        case 0 % existing
            ncm.mu = value(mu);
            % ncm.X = value(X);            
            ncm.X = cond(ncm.M);

            ncm.W_bar = value(W_bar);
            ncm.optDone = sol.problem;

            ncm.M = inv(ncm.W_bar/ncm.mu);

        case 1 % proposed 1 (effective space)
            ncm.M = value(M);
            ncm.mu = value(mu);
            ncm.X = cond(ncm.M);
            ncm.optDone = sol.problem;

        case 2 % proposed 2 (inverse)
            % ncm.AUX = value(AUX);
            % ncm.mu = value(mu);
    end

    %% CONTRACTION CONDITION CHECK
    switch ncm.ctrl_no 
        case 0
            M = inv(ncm.W_bar/ncm.mu); pre_M = inv(pre_W_bar/pre_mu);
            cond_check = (M-pre_M)/dt + (SDC'*M + M*SDC) - 2*M*B*inv_R*B'*M + 2*alpha*M;
        case 1 % proposed 1 (effective space)
            cond_check = (ncm.M-pre_M)/dt + (SDC'*ncm.M + ncm.M*SDC) - 2*ncm.M*B*inv_R*B'*ncm.M + 2*alpha*ncm.M;
        case 2 % proposed 2 (inverse)
            error("not maintained anymore");
    end
    ncm.contraction_flag = any(eig(cond_check) <= 0);
end