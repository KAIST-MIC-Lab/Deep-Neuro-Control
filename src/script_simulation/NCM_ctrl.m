function ncm = NCM_ctrl(ncm, x, xd, ud, param)
    
    dt = ncm.dt;
    d_MAX = ncm.d_MAX;  % maximum disturbance
    m_MIN = ncm.m_MIN;

    mu = ncm.mu; 

    alpha = ncm.alpha;  % dacay rate

    inv_R = ncm.inv_R;  

    x_num = ncm.x_num;

    %% 
    th = x(1);
    L = param.L;    % Inductance (mH)
    P = param.P;         % Pole pairs
    Phi = param.Phi; % Flux (Wb)
    B = [0 0; -(3*P*Phi)/(2*L)*sin(P*th) (3*P*Phi)/(2*L)*cos(P*th)]; 

    SDC = SDCmotor(x, xd, ud, param);  % System Dynamics Coefficient

    %% OPTIMIZATION PROBLEM DEFINITION
    W = sdpvar(x_num ,x_num);
    assign(W, eye(x_num)*0);
    X = sdpvar(1,1);
    % mu = sdpvar(1,1);
    assign(X, 1);
    % assign(mu, 1);

    c0 = 1;
    c1 = 0;

    pre_W = ncm.W;

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

        ncm.init = 0;
    end

    ops = sdpsettings('solver', 'sedumi', ...
                  'verbose', 0, ...
                  'sedumi.eps', 1e-8);
    % optimize(con, obj);
    sol = optimize(con, obj, ops);

    if sol.problem ~= 0
        ncm.optDone = 0;
        ncm.X = 0;
        ncm.W = pre_W;  % keep the previous value
        % warning('YALMIP Error: %s', yalmiperror);
    else
        ncm.W = value(W);
        ncm.X = value(X);
        ncm.optDone = 1;
    end
    
    %%

end