    function ncm = NCM_ctrl(ncm, x, xd, ud, param)
    % for initial control
    if ncm.init || ncm.ctrl_no == 0 
        ncm.optDone = 0;
        ncm.init = 0;
        ncm.contraction_flag = 1;
        ncm.u = ud;
        ncm.cmp_t = 0;
        return
    end

    %% FROM STRUCTURE TO LOCAL VARIABLES
    dt = ncm.dt;
    x_num = ncm.x_num;
    e = x-xd;
    u_max = param.max_u;

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
        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 1 
            % prepare
            pre_M = ncm.M; 
            pre_y = ncm.y;
            pre_s = ncm.nu;
            M_energy = ncm.M_energy;
            
            % pre-calculation
            H = (-2*B*inv_R*B');
            f = transpose(e'*(1/dt*eye(x_num)+2 *SDC'+2*alpha*eye(x_num)));
            c = -(e'*pre_M*e)/dt;

            Hc = (B*inv_R*inv_R*B');
            fc = transpose(-2*ud'*inv_R*B');
            cc = (ud'*ud - u_max^2);

            max_Ctr_y = -1/2 * H\f;
            max_Ctr = max_Ctr_y'*H*max_Ctr_y + f'*max_Ctr_y + c; 

            e_norm = norm(e,2);
            e_norm_sq = e'*e;

            % optimization variable
            y = sdpvar(x_num,1); assign(y, pre_y);
            s = sdpvar(1,1); assign(s, pre_s);

            % objective function
            obj = (e_norm*norm(y,2) - e'*y) + ncm.lbd*s;

            % constraint
            con = [
                (e'*y - M_energy*e_norm_sq)^2 <= (1e-2+s)*e_norm_sq^2; % energy constraint
                s >= 0; % slack variable
                y'*H*y + f'*y + c <= 0; % contraction constraint
                (y'*Hc*y + fc'*y + cc)/abs(1e2) <= 0; % saturation constraint
                % e'*y >= 1e-5;
            ];

        % EXISTING – CV-STEM
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 3 
            % prepare
            pre_W_bar = ncm.W_bar; 
            pre_nu = ncm.nu;

            % optimization variable 
            % (for SEMUDI solver, warm-start not supported)
            W_bar = sdpvar(x_num ,x_num); % assign(W_bar, pre_W_bar);
            nu = sdpvar(1,1); % assign(nu, pre_nu);
            X = sdpvar(1,1); % assign(X, ncm.X);

            % objective function
            obj = X + ncm.lbd * nu;
        
            % constraint
            con = [
                -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*nu*B*inv_R*B' <= -2*alpha*W_bar;
                eye(x_num) <= W_bar;
                W_bar <= X * eye(x_num);
                % X == nu/m_min;
                % m_min <= nu;
                1e-99 <= nu;
                1 <= X;
            ];
    end
    
    %% OPTIMIZATION 
    switch ncm.ctrl_no
        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 1
            ops = sdpsettings('solver', 'fmincon');
            ops = sdpsettings(ops, 'fmincon.Algorithm', 'sqp');
            ops = sdpsettings(ops, 'verbose', 0);
            ops = sdpsettings(ops, 'debug',0);
            ops = sdpsettings(ops, 'usex0', 1);

            ops = sdpsettings(ops, 'fmincon.MaxIter', 1e9);
            % ops = sdpsettings(ops, 'fmincon.TolFun', 1e-5); % obj tol
            % ops = sdpsettings(ops, 'fmincon.TolCon', 1e-5); % cstr tol 
            % ops = sdpsettings(ops, 'fmincon.TolX', 1e-7); % step tol

        % EXISTING – CV-STEM
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 3
            ops = sdpsettings('solver', 'sedumi');
            ops = sdpsettings(ops, 'verbose', 0);
            ops = sdpsettings(ops, 'debug',0);
            % ops = sdpsettings(ops, 'usex0', 1);
    
    end

    % optimize!
    tic;
    sol = optimize(con, obj, ops);
    ncm.cmp_t = toc;

    % optimization result check 
    %   (0: success, 1: infeasible, 2: unbounded, 3: max iterations, 4: numerical error)
    %   see https://yalmip.github.io/command/yalmiperror/

    %% TERMINIATING OPTIMIZATION
    ncm.optDone = sol.problem;

    switch ncm.ctrl_no
        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 1 
            if sol.problem == 0 || norm(e) > 1e-9 
                if value(y)'*e < -1e-5
                    error("y'*e is negative")
                end

                ncm.y = value(y);
                ncm.nu = value(s); % <======== CORRECT!!!

                % condition number and metric
                th = acos(e'*ncm.y/(norm(e)*norm(ncm.y)));
                th = real(th); % numerical correction
                X = (1+sin(th))/(1-sin(th));

                M = diag([sqrt(X);ones(x_num-2,1);1/sqrt(X)]);
                ncm.M = (M/norm(M)) * (norm(ncm.y)/norm(e));

                if ~isreal(ncm.M)
                    fprintf("check")
                end

                ncm.X = cond(ncm.M);

                % BFGS
                % ncm.M = eye(x_num) ... 
                %     - (e*e')/(e'*e) ...
                %     + (ncm.y*ncm.y')/(ncm.y'*e);
                % ncm.M = pre_M ... 
                %     - ((pre_M*e)*(e'*pre_M))/(e'*pre_M*e) ...
                %     + (ncm.y*ncm.y')/(ncm.y'*e);

            end
            ncm.u = ud - ncm.inv_R*B' * ncm.y;

        % EXISTING – CV-STEM
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 3
            if sol.problem == 0
                ncm.nu = value(nu);
                ncm.W_bar = value(W_bar);
                ncm.M = ncm.nu * (ncm.W_bar\eye(size(ncm.W_bar)));  % M = nu * inv(W_bar)
                ncm.X = cond(ncm.M);
            end
            ncm.u = ud - ncm.inv_R*B' * ncm.M * e;
    end

    %% CONTRACTION CONDITION CHECK
    switch ncm.ctrl_no 
        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 1 
            con_check_vector = ncm.y'*H*ncm.y + f'*ncm.y + c;
            cond_check = (ncm.M-pre_M)/dt + (SDC'*ncm.M + ncm.M*SDC) - 2*ncm.M*B*inv_R*B'*ncm.M + 2*alpha*ncm.M;
            max_eig = max(eig(cond_check));
            % 
            % fprintf("c_vec: %d, c: %d, c*e: %d\n", con_check_vector<=1e-5, max_eig<=1e-5, e'*cond_check*e <= 1e-5);
            % 
            % if max_eig > 1e-9
            %     fprintf("chec")
            % end

            ncm.contraction_flag = con_check_vector <= 1e-9;

        % EXISTING – CV-STEM
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 3 
            pre_M = pre_nu*pre_W_bar\eye(size(pre_W_bar));
            cond_check = (ncm.M-pre_M)/dt + (SDC'*ncm.M + ncm.M*SDC) - 2*ncm.M*B*inv_R*B'*ncm.M + 2*alpha*ncm.M;
            max_eig = max(eig(cond_check));

            ncm.contraction_flag = max_eig <= 1e-9;

    end
    

end
