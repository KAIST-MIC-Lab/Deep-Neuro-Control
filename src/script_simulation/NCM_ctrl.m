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
        case 3 % existing
            % prepare
            pre_W_bar = ncm.W_bar; 
            pre_mu = ncm.mu;
            e = x-xd;

            % optimization variable
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


        case 12 % proposed 1 (effective space)
            % prepare
            e = x-xd;
            pre_M = ncm.M; 
            pre_y = ncm.y;
            M_energy = ncm.M_energy;

            % optimization variable
            y = sdpvar(x_num,1);
            eps = sdpvar(1,1); assign(eps, 0);

            % optimal point (with respect to 1 cond)
            H = (-2*B*inv_R*B');
            f = transpose(e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num)));
            c = -(e'*pre_M*e)/dt;
            
            % mu_h = e'*H*e; mu_f = f'*e; % c = 0;
            % mu_trg_1 = 1/2/mu_h * (-mu_f-sqrt(mu_f^2-4*mu_h*c));
            % % mu_trg_1 = 1/2/mu_h * (-mu_f+sqrt(mu_f^2-4*mu_h*c));
            % if mu_f^2-4*mu_h*c < 0
            % %     y_trg = zeros(x_num,1);
            %     warning("imaginary root")
            % end
            %     % else
            % mu_trg_1 = real(mu_trg_1);
            % y_trg = mu_trg_1 * e;
            % % end

            assign(y, pre_y);
            
            % objective function
            lbd = 1e-1; % slack var penalty
            obj = (sqrt((e'*e)*(y'*y)) - e'*y) + lbd*eps;
            % obj = (y-y_trg)'*(y-y_trg);
            % obj = (ud - inv_R*B'*y)'*(ud - inv_R*B'*y);

            % constraint
            u_max = param.max_u;
            con = [
                (e'*y - M_energy*(e'*e))^2 <= (1e-5+eps)*(e'*e);
                eps >= 0
                % y'*H*y + f'*y <= 0;
                y'*H*y + f'*y + c <= 0;
                y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
                % y'*e >= 1e-9;
            ];
    end
    
    %% OPTIMIZATION 
    switch ncm.ctrl_no
        case 3
            ops = sdpsettings('verbose', 0);
            ops = sdpsettings(ops, 'debug',0);
    
            % ops = sdpsettings(ops, 'solver','mosek');
            ops = sdpsettings(ops,'solver', 'sedumi');
            % ops = sdpsettings(ops,'sedumi.eps', 1e-1);
            % ops = sdpsettings(ops, 'sedumi.cg.maxiter', 1000);
            
        case {13, 14, 12}
            ops = sdpsettings('solver', 'fmincon', ...
                'verbose', 0, ...
                'fmincon.Algorithm', 'sqp', ...
                'fmincon.TolFun', 1e-8, ... % obj func 
                'fmincon.MaxIter', 1000000000);
                % 'MaxIterations', 500000);
            ops = sdpsettings(ops, 'fmincon.TolCon', 1e-8);
            ops = sdpsettings(ops, 'debug',0);

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
        case 3 % existing
            if sol.problem == 0
                ncm.mu = value(mu);
                ncm.W_bar = value(W_bar);
                ncm.M = ncm.mu * (ncm.W_bar\eye(size(ncm.W_bar)));  % M = mu * inv(W_bar)
                ncm.X = cond(ncm.M);

                ncm.optDone = sol.problem;
    
            end

            ncm.u = ud - ncm.inv_R*B' * ncm.M * e;

        case 12 % proposed 1 (effective space)

            if sol.problem == 0
                ncm.y = value(y);
                ncm.mu = eps; % <======== CORRECT!!!

                % pre_res = check(con);
                % if pre_res(end) <= 0
                %     ncm.M = mu_trg_1*eye(x_num);
                % else
                    % if norm(ncm.y) == 0 || norm(e) == 0
                    %     tmp1 = 0;
                    %     tmp2 = 0;
                    % else
                        tmp1 = - ((pre_M*e)*(e'*pre_M))/(e'*pre_M*e);
                        tmp2 = + (ncm.y*ncm.y')/(ncm.y'*e);
                    % % end
                    ncm.M = pre_M + tmp1 + tmp2;

                % ncm.M = eye(x_num) - (e*e')/(e'*e) + (ncm.y*ncm.y')/(ncm.y'*e);
                % end
                ncm.X = cond(ncm.M);
            end
            
            % cond0 = ncm.y'*H*ncm.y+f'*ncm.y+c <= 0;
            % fprintf("y: %.3f %.3f %.3f, cond0: %d\n", ncm.y(1), ncm.y(2), ncm.y(3), cond0);
            % fprintf("eps: %.4f\n", value(eps));

            ncm.u = ud - ncm.inv_R*B' * ncm.y;

    end

    %% CONTRACTION CONDITION CHECK
    switch ncm.ctrl_no 
        case 3 % existing
            pre_M = pre_mu*pre_W_bar\eye(size(pre_W_bar));
            cond_check = (ncm.M-pre_M)/dt + (SDC'*ncm.M + ncm.M*SDC) - 2*ncm.M*B*inv_R*B'*ncm.M + 2*alpha*ncm.M;
        case {13, 14, 12} % proposed 1 (effective space)
            cond_check = (ncm.M-pre_M)/dt + (SDC'*ncm.M + ncm.M*SDC) - 2*ncm.M*B*inv_R*B'*ncm.M + 2*alpha*ncm.M;

    end
    ncm.contraction_flag = all(eig(cond_check) <= 0);

    %%
    return

end

        % case 12 % proposed 1 (effective space)
        %     pre_M = ncm.M; 
            
        %     pre_y = ncm.y;
        %     e = x-xd;

        %     M = sdpvar(x_num, x_num); assign(M, pre_M);
        %     mu = sdpvar(1,1); assign(mu, ncm.mu);
        %     y = sdpvar(x_num,1); assign(y, pre_y);


        %     % objective function
        %     obj = norm(M-mu*eye(x_num)) - 12e-1*mu;

        %     % constraint
        %     u_max = param.max_u;
        %     con = [
        %         -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num))*y <= 0;
        %         % -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num))*y + (e'*pre_M*e)/dt <= 0;
        %         % y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
        %         M*e == y;
        %         mu >= 1e-9;
        %     ];

        % case 11 % [BACKUP] proposed 1 (effective space)
        %     pre_M = ncm.M; 
        %     % pre_M = zeros(size(ncm.M)); 
        %     % pre_M_bar = ncm.M_bar;
        %     pre_X = ncm.X;
        %     pre_y = ncm.y;
        %     e = x-xd;

        %     M = sdpvar(x_num, x_num); assign(M, pre_M);
        %     % M_bar = sdpvar(x_num, x_num); assign(M_bar, pre_M_bar);
        %     mu = sdpvar(1,1); assign(mu, ncm.mu);
        %     % X = sdpvar(1,1); assign(X, pre_X);

        %     y = sdpvar(x_num,1); assign(y, pre_y);

        %     % mu_min = sdpvar(1,1);
        %     % mu_max = sdpvar(1,1);

        %     % objective function
        %     % obj = norm(M-eye(x_num));
        %     obj = norm(M-mu*eye(x_num)) - 1e0*mu;
        %     % obj = norm(M_bar-eye(x_num)) - mu;
        %     % obj = norm(M - (trace(M)/x_num)*eye(x_num)) - (trace(M)/x_num) * 1e5;
        %     % obj = norm(M - (trace(M)/x_num)*eye(x_num), 'fro')^2 - (trace(M)/x_num) * 1e1;
        %     % obj = log(trace(M)) - 1/x_num * log(det(M));
        %     % obj = -mu;
        %     % obj = -mu^2;
        %     % obj = 1/X;

        %     m_min = 1e-99;
        %     m_max = 1e9;

        %     % constraint
        %     u_max = param.max_u;
        %     con = [
        %         -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num))*y - (e'*pre_M*e)/dt <= 0;
        %         y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
        %         M*e == y;
        %         % M >= mu*eye(x_num);
        %         % M <= m_max * eye(x_num);
        %         % M >= m_min * eye(x_num);
        %         % M <= X*mu * eye(x_num);
        %         mu >= 1e-9;
        %         % M-pre_M <= 1e-3*eye(x_num);
        %         % -M+pre_M <= 1e-3*eye(x_num);
        %         % X >= 1e-99;
        %         % M == mu * M_bar;
        %         % e'*y >= e'*pre_M*e;
        %         % M >= 1e-99 * eye(x_num);
        %         % M <= 1e2*eye(x_num);
        %         % X * 1e2 * eye(x_num) <= M;
        %         % X >= 1e-1;
        %         % X <= 1;
        %     ];

        % case 2 % proposed 2 (inverse)
        %     % AUX = sdpvar(x_num, x_num);
    
        %     % mu = sdpvar(1,1);
        %     % assign(mu, ncm.mu);
    
        %     % assign(AUX, ncm.AUX);

        %     % % objective function
        %     % lbd = 1e-2;
        %     % % obj = X + trace(AUX);
        %     % obj = X + trace(AUX)+ lbd * mu;

        %     % % constraints
        %     % con = [
        %     %     -(W_bar-pre_W_bar)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
        %     %     eye(x_num) <= W_bar;
        %     %     W_bar <= X * eye(x_num);
        %     %     % X == mu/m_min;
        %     %     % m_min <= mu;
        %     %     1e-99 <= mu;
        %     %     1 <= X;
        %     % ];
        
        %     % max_u = param.max_u;
        %     % con = [
        %     %     [max_u^2  (W_bar*inv(inv_R*B')*ud-mu*(x-xd))'
        %     %     (W_bar*inv(inv_R*B')*ud-mu*(x-xd)) AUX] >= 0;
        %     %     [AUX (W_bar*inv(inv_R*B'))'
        %     %     (W_bar*inv(inv_R*B')) eye(x_num)] >= 0,
        %     %     con, 
        %     % ];

% ======================================================================

        % case 14 % proposed 1 (effective space)
        %     pre_mu = ncm.mu;
        %     pre_M = ncm.M - ncm.mu*eye(x_num);
        %     e = x-xd;
        %     for idx = 1:length(e)
        %         if abs(e(idx)) < 1e-9
        %             e(idx) = sign(e(idx)) * 1e-9;
        %         end
        %     end

        %     mu = sdpvar(1,1); assign(mu, 1e4);
        %     del_M = sdpvar(x_num, x_num); assign(del_M, zeros(x_num,x_num))
        %     del_y = sdpvar(x_num,1); assign(del_y, zeros(x_num,1))
        %     % del_y = zeros(3,1);
        %     % eps = sdpvar(1,1);
        %     eps = zeros(1,1);

        %     % constraint
        %     H = [
        %         -2*e'*B*inv_R*B'*e, -2*e'*B*inv_R*B';
        %         -2*B*inv_R*B'*e, -2*B*inv_R*B'
        %     ];
        %     f = [
        %         e'*e/dt+e'*SDC'*e+e'*SDC*e+2*alpha*(e'*e);
        %         transpose(e'/dt+2*e'*SDC'+2*alpha*e')    
        %     ];
        %     c = -(e'*e)*pre_mu/dt + -e'*del_M*e/dt;

        %     % objective function
        %     % obj = norm(del_y);
        %     obj = del_y'*del_y + 5e-1*eps;
        %     % obj = del_y'*del_y + eps -1e2*mu;
        %     % lbd = 0e12;
        %     % obj = lbd*(del_y'*del_y) + ([mu;del_y]'*H*[mu;del_y] + f'*[mu;del_y])'*([mu;del_y]'*H*[mu;del_y] + f'*[mu;del_y]); % -> a m + b <= 0
        %     % obj =  lbd*(del_y'*del_y) +([mu;del_y]- .5*inv(H)*f)'*([mu;del_y]- .5*inv(H)*f);
        %     % obj =  lbd*(del_y'*del_y) +([mu;del_y]- H\f)'*([mu;del_y]- H\f);

        %     u_max = param.max_u;
        %     con = [
        %         % -2*(e'*B*inv_R*B'*e)*mu^2 -4*del_y'*(B*inv_R*B')*mu*e -2*del_y'*(B*inv_R*B')*del_y ...
        %         % + (e'/dt+e'*SDC'+e'*SDC+2*alpha*e')*(mu*e+del_y) <= -0;
        %         % mu >= mu_bound;
        %         % -2*e'*B*inv_R*B'*e*mu^2 + (e'*e/dt+e'*SDC'*e+e'*SDC*e+2*alpha*(e'*e))*mu <= 0;
        %         [mu;del_y]'*H*[mu;del_y] + f'*[mu;del_y] + - eps <= 0; % -> a m + b <= 0
        %         % [mu;del_y]'*H*[mu;del_y] + f'*[mu;del_y] + c <= 0;
        %         % H*[mu;del_y]]
        %         % mu^2*a + mu*(b+2*alpha*(e'*e)) +c <= 0
        %         % e'*(mu*eye(x_num)+del_M)*(B*inv_R*inv_R*B')*(mu*eye(x_num)+del_M)*e + (-2*ud'*inv_R*B')*(mu*eye(x_num)+del_M)*e + (ud'*ud - u_max^2) <= 0;
        %         (e'*mu+del_y')*(B*inv_R*inv_R*B')*(mu*e+del_y) + (-2*ud'*inv_R*B')*(mu*e+del_y) + (ud'*ud - u_max^2) <= 0;
        %         % del_y == del_M * e;
        %         % mu >= 1e-9;
        %         mu >= 2.5e1
        %         del_y'*e >= 0
        %         del_y'*e <= 0.0001
        %         eps >= 0
        %         % del_M >= 0;
        %         % mu <= 20;
        %     ];

        % case 13 % proposed 1 (effective space)
        %     pre_mu = ncm.mu;
        %     e = x-xd;
        %     for idx = 1:length(e)
        %         if abs(e(idx)) < 1e-9
        %             e(idx) = sign(e(idx)) * 1e-9;
        %         end
        %     end

        %     mu = sdpvar(1,1); assign(mu, ncm.mu);
            
            
        %     % constraint
        %     a = -2*e'*B*inv_R*B'*e;
        %     b = 1/dt*(e'*e) + e'*SDC'*e + e'*SDC*e;
        %     c = -(e'*e)*pre_mu/dt;
            
        %     mu_bound = 1/a * (-b-2*alpha*(e'*e));
        %     fprintf("mu bound %.4f\n", mu_bound);
            
        %     % objective function
        %     obj = (mu-mu_bound)^2; % mu > mu_bound, mu, minimum mu

        %     u_max = param.max_u;
        %     con = [
        %         % mu >= mu_bound;
        %         % mu^2*a + mu*(b+2*alpha*(e'*e)) +c <= 0
        %         mu^2*(e'*(B*inv_R*inv_R*B')*e) + (-2*ud'*inv_R*B')*e*mu + (ud'*ud - u_max^2) <= 0;
        %         mu >= 1e-9;
        %         % mu <= 20;
        %     ];