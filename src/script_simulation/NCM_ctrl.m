function ncm = NCM_ctrl(ncm, x, xd, ud)
    
    dt = ncm.dt;
    d_MAX = ncm.d_MAX;  % maximum disturbance
    
    mu = ncm.mu; 

    alpha = ncm.alpha;  % dacay rate

    inv_R = ncm.inv_R;  

    %% 
    [f, B] = system_func(x);
    [fd, Bd] = system_func(xd);

    SDC = getSDC(x, xd);  % System Dynamics Coefficient

    % try
    %     % SDC = (f+B*ud - (fd+Bd*ud)) \ (x - xd);
    %     for i = 1:4
    %         for j = 1:4
    %             SDC(i,j) = (f(i) + B(i,:)*ud - (fd(j) + Bd(j,:)*ud)) / (x(i) - xd(j));
    %             if isinf(SDC(i,j))
    %                 SDC(i,j) = 0.0001;
    %             end
    %         end
    %     end
    % catch
    %     error('System is not controllable');
    % end

    %% OPTIMIZATION PROBLEM DEFINITION
    W = sdpvar(4,4);
    assign(W, eye(4)*1);
    X = sdpvar(1,1);
    assign(X, 1);

    pre_W = ncm.W;

    obj = X;
    % if sum(pre_W, 'all') == 0
        % con = [
        %     % W == W';
        %     mu*(W'*SDC' + SDC*W) - mu*2*d_MAX*B*inv_R*B' <= -mu*2*alpha*W;
        %     eye(4) <= mu*W;
        %     mu*W <= X * eye(4);
        % ];

        % con = [
        %     [
        %     mu*(W'*SDC' + SDC*W)-mu*2*d_MAX*B*inv_R*B'+mu*2*alpha*W mu*W;
        %     mu*W zeros(4,4)
        %     ] <= 0
        %     eye(4) <= mu*W;
        %     mu*W <= X * eye(4);    
        % ];
    % else
        con = [
            % W == W';
            -mu*(W-pre_W)/dt + mu*(W'*SDC' + SDC*W) - mu*2*d_MAX*B*inv_R*B' <= -mu*2*alpha*W;
            eye(4) <= mu*W;
            mu*W <= X * eye(4);
        ];
    % end

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

function [f, B] = system_func(x)
    th = x(1);  % theta
    thd = x(2);  % theta dot
    ia = x(3);  % Ia
    ib = x(4);  % Ib
    
    %%
    L = .66e-3;    % Inductance (mH)
    R = 0.251;     % Resistance (Ohm)
    J = 3.24e-5;   % Inertia (kg.m^2)
    Phi = 16.8e-3; % Flux (Wb)
    P = 4;         % Pole pairs
    fv = 2e-3;     % Viscous friction (N.m.s/rad)
    
    %% 
    trq = -(3/2)*P*Phi*sin(P*th)*ia + (3/2)*P*Phi*cos(P*th)*ib;

    %%
    f = [
        thd;
        (1/J)*(trq - fv*thd);
        (1/L)*(- R*ia + P*Phi*thd*sin(P*th));
        (1/L)*(- R*ib - P*Phi*thd*cos(P*th))
    ];

    B = [
        0, 0;
        0, 0;
        1/L, 0;
        0, 1/L
    ];

end

function SDC = getSDC(x,xd)

    %%
    L = .66e-3;    % Inductance (mH)
    R = 0.251;     % Resistance (Ohm)
    J = 3.24e-5;   % Inertia (kg.m^2)
    Phi = 16.8e-3; % Flux (Wb)
    P = 4;         % Pole pairs
    fv = 2e-3;     % Viscous friction (N.m.s/rad)

    %%
    a21 = 1/J * (-(3/2)*P*Phi) * ...
        (-intXcosY(x(3),xd(3), x(1),xd(1), P) + ...
        intXsinY(x(4),xd(4), x(1),xd(1), P));
    a23 = 1/J * (-(3/2)*P*Phi) * ...
        intXsinY(1,1, x(1),xd(1), P);
    a24 = 1/J * (+(3/2)*P*Phi) * ...
        intXcosY(1,1, x(1),xd(1), P);

    a31 = 1/L*P^2*Phi * ...
        intXcosY(x(2),xd(2), x(1),xd(1), P);
    a41 = 1/L*P^2*Phi * ...
        intXsinY(x(2),xd(2), x(1),xd(1), P);

    a32 = 1/L*P*Phi * ...
        intXsinY(1,1, x(1),xd(1), P);
    a42 = 1/L*P*Phi * ...
        intXcosY(1,1, x(1),xd(1), P);

    %%
    SDC = [
        0 1 0 0;
        a21 -fv/J a23 a24;
        a31 a32 0 -R/L;
        a41 a42 -R/L 0
    ];

    for i = 1:1:4
        for j = 1:1:4
            if isnan(SDC(i,j))
                SDC(i,j) = 0;
            end
        end
    end
end

function z = intXsinY(x1,x2, y1,y2, a)
    z = 1/(a*(y1-y2)) * (-cos(a*y1)*x1 + cos(a*y2)*x2) + (x1-x2)/(a*(y1-y2))^2 * (sin(a*y1) - sin(a*y2));
end

function z = intXcosY(x1,x2, y1,y2, a)
    z = 1/(a*(y1-y2)) * (sin(a*y1)*x1 - sin(a*y2)*x2) + (x1-x2)/(a*(y1-y2))^2 * (-cos(a*y1) + cos(a*y2));
end