function ncm = NCM_ctrl2(ncm, x, xd, ud, param)
% NCM_ctrl2 (YALMIP + SDPT3, robustified)
% - Slack s >= 0 for LMI relaxation with heavy penalty
% - Warm-start + symmetric/PSD post-fix
% - Bounds guards for W_bar/mu/X
% - One adaptive retry with larger eps if solver fails

    % (선택) YALMIP 캐시 클리어: 문제에 따라 켜/끌기
    % yalmip('clear');

    if ncm.init
        ncm.optDone = 0;
        ncm.init = 0;
        return
    end

    %% Unpack & defaults
    dt     = ncm.dt;
    n      = ncm.x_num;
    alpha  = ncm.alpha;
    inv_R  = ncm.inv_R;
    lbd    = ncm.lbd;

    if ~isfield(ncm,'mu_min'),    ncm.mu_min    = 1e-4; end
    if ~isfield(ncm,'mu_max'),    ncm.mu_max    = inf;  end
    if ~isfield(ncm,'eps_reg'),   ncm.eps_reg   = 1e-8; end       % metric reg (for control law)
    if ~isfield(ncm,'X_min'),     ncm.X_min     = 1.0;  end
    if ~isfield(ncm,'X_max'),     ncm.X_max     = inf;  end
    if ~isfield(ncm,'eps_psd0'),  ncm.eps_psd0  = 1e-8; end       % LMI PSD margin (first try)
    if ~isfield(ncm,'eps_psd1'),  ncm.eps_psd1  = 1e-6; end       % LMI PSD margin (retry)
    if ~isfield(ncm,'rho_slack'), ncm.rho_slack = 1e3;  end       % slack penalty

    %% System matrices
    th = x(1);
    J  = param.J;  P = param.P;  Phi = param.Phi;

    B = [ 0 0;
         -(3*P*Phi)/(2*J)*sin(P*th),  (3*P*Phi)/(2*J)*cos(P*th) ];

    SDC = SDCmotor(x, xd, ud, param);

    % Finite check
    if any(~isfinite([x(:); xd(:); ud(:); B(:); SDC(:)]))
        warning('NCM_ctrl2: Non-finite detected (inputs/B/SDC). Keeping previous values.');
        ncm.optDone = -99;
        return
    end

    %% Decision variables
    W_bar = sdpvar(n,n,'symmetric');
    X     = sdpvar(1,1);
    mu    = sdpvar(1,1);
    s     = sdpvar(1,1);      % slack ≥ 0

    aux   = sdpvar(2,1);
    % W     = sdpvar(n,n,'symmetric');

    % Warm start
    pre_W = ncm.W_bar;
    % assign(W_bar, pre_W);
    % assign(X,     max(ncm.X, ncm.X_min));
    % assign(mu,    max(ncm.mu, ncm.mu_min));
    % assign(s,     0);

    %% One-shot builder for constraints & objective (eps_psd param)
    function [Con,Obj] = build_problem(eps_psd)
        Con = [];
        % Main LMI with slack: LHS <= RHS + s*I
        %   LHS = -(W - pre_W)/dt + (W*SDC' + SDC*W) - 2*mu*B*inv_R*B'
        %   RHS = -2*alpha*W - eps_psd*I
        LHS = -(W_bar - eye(n))/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*(B*inv_R*B');
        RHS = -2*alpha*W_bar - eps_psd*eye(n);
        Con = [Con,  LHS <= RHS + s*eye(n) ];

        % Control input norm constaint
        % u_bar = 0.1; e = x-xd;
        % 
        % cstr_Mat = [
        %     2*u_bar^2-ud'*ud--2*(ud'*inv_R*B')*aux     aux'
        %     aux                                         inv(B*inv_R*inv_R*B')
        % ];

        % Con = [Con, 0 <= cstr_Mat, e==W_bar*aux];

        % Con = [Con, W_bar*inv(inv_R*B')*(ud-u_bar)<=e/mu];

        % Bounds and PSD guards
        Con = [Con,  W_bar >= ncm.X_min*eye(n) ];
        if isfinite(ncm.X_max), Con = [Con, W_bar <= X*eye(n) ]; else, Con = [Con, W_bar <= X*eye(n) ]; end
        Con = [Con,  X >= ncm.X_min ];
        if isfinite(ncm.X_max), Con = [Con, X <= ncm.X_max ]; end

        Con = [Con,  mu >= ncm.mu_min ];
        if isfinite(ncm.mu_max), Con = [Con, mu <= ncm.mu_max ]; end

        Con = [Con,  s >= 0 ];

        % Objective: minimize X + lbd*mu + rho*s (slack heavily penalized)
        Obj = X + lbd*mu + ncm.rho_slack*s;
    end

    %% Solve helper
    function ok = solve_with_eps(eps_psd)
        [Con,Obj] = build_problem(eps_psd);
        opts = sdpsettings('solver','sdpt3', ...
                           'verbose',0,'debug',1, ...
                           'savesolverinput',1,'savesolveroutput',1, ...
                           'cachesolvers',0);
        try
            sol = optimize(Con, Obj, opts);
        catch ME
            warning('SDPT3 crashed: %s', ME.message);
            ncm.optDone = 999;
            ok = false; return
        end
        if sol.problem ~= 0
            warning('YALMIP/SDPT3 failed (eps=%.1e): %s', eps_psd, sol.info);
            ncm.optDone = sol.problem;
            ok = false; return
        end
        % Update (store in caller)
        ncm.X  = max(value(X),  ncm.X_min);
        ncm.mu = max(value(mu), ncm.mu_min);
        Wn = value(W_bar);
        Wn = 0.5*(Wn + Wn');               % sym fix

        % PSD post-fix
        [V,D] = eig((Wn+Wn')/2);
        D = max(D, eps_psd*eye(n));
        Wn = V*D*V';

        % Clamp upper bound if X_max < inf
        if isfinite(ncm.X_max)
            % ensure Wn <= X*I by clipping eigenvalues
            lam = diag(D);
            lam = min(lam, ncm.X_max);
            Wn = V*diag(lam)*V';
        end

        if any(~isfinite(Wn(:))) || ~isfinite(ncm.mu) || ~isfinite(ncm.X)
            warning('NCM_ctrl2: Non-finite solution. Keeping previous values.');
            ncm.optDone = -97;
            ok = false; return
        end

        ncm.W_bar  = Wn;
        ncm.optDone = 0;
        ok = true;
    end

    %% Try solve with eps_psd0, then eps_psd1 if needed
    if ~solve_with_eps(ncm.eps_psd0)
        % Retry once with looser PSD margin
        solve_with_eps(ncm.eps_psd1);
    end
end
