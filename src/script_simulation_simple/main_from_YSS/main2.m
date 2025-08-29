clearvars -global
clear all %#ok<CLALL>
close all

% YALMIP/solvers path
%% PATH
addpath(genpath([pwd filesep 'YALMIP-master']));
% addpath(genpath([pwd filesep 'sdpt3-master']));
addpath(genpath([pwd filesep 'SDPT3-4.0']));
%% FASTEN YOUR SEATBELT
clear
RESULT_SAVE_FLAG = 0;   % save .mat
FIGURE_PLOT_FLAG = 1;   % plot
FIGURE_SAVE_FLAG = 0;   % save figs

%% SIMULATION SETTING
T       = 2;          % [s]
ctrl_dt = 1e-3;       % [s] (1 kHz)
dt      = 1e-4;       % [s] (1 MHz)  -> 매우 촘촘: 필요시 1e-5~1e-4 권장
rpt_dt  = 1e-1;       % [s]
t       = 0:dt:T;

% controller/report stride
ctrl_stride = max(1, round(ctrl_dt/dt));
rpt_stride  = max(1, round(rpt_dt/dt));

I_MAX  = 2;    % ★ 0.5A로 강하게 제한 (필요시 0.3~0.8A)
du_max = 0.05;   % ★ 1 step당 0.01A만 변하게 (slew-rate ~10 A/s @1kHz)
u_prev = [0;0];
ncm.kc=1;
% (옵션) ud 필터
ud_f     = [0;0];
tau_ud   = 0.5;                       % 50 ms
alpha_ud = min(1, ctrl_dt/tau_ud);
y_lp   = 0;
%% SM PARAMETERS
param.L   = .66e-3;     % H (unused here)
param.R   = 0.251;      % Ohm (unused here)
param.J   = 3.24e-5;    % kg.m^2
param.Phi = 16.8e-3;    % Wb
param.P   = 4;          % -
param.fv  = 2e-3;       % N.m.s/rad

%% REPORT
fprintf("\n*** SIMULATION INFORMATION ***\n");
fprintf("Termination Time : %.2f s\n", T);
fprintf("Controller dt    : %.2e s\n", ctrl_dt);
fprintf("Simulation dt    : %.2e s\n", dt);
fprintf("Report dt        : %.2e s\n\n", rpt_dt);
fprintf("RESULT_SAVE_FLAG : %d\n", RESULT_SAVE_FLAG);
fprintf("FIGURE_PLOT_FLAG : %d\n", FIGURE_PLOT_FLAG);
fprintf("FIGURE_SAVE_FLAG : %d\n\n", FIGURE_SAVE_FLAG);

%% STATES / INPUTS / REFS
x     = [0;0];
x_non = [0.3;0];
u     = [0;0];
xd    = [0;0];

% reference (theta only)
xd_func   = @(tt) 2*(cos(tt)-1) + 1*(cos(10*tt)-1)+0.3;
dxd_func  = @(tt) -2*sin(tt)    - 10*sin(10*tt);
ddxd_func = @(tt) -2*cos(tt)    - 100*cos(10*tt);

% dimensions
num_x = numel(x);
num_u = numel(u);
num_t = numel(t);

%% CONTROLLER INIT (stabilized defaults)
NCM_init     % ← 이전에 제공한 안정화 프리셋 사용

if ~isfield(ncm,'mu_min'),  ncm.mu_min = 1e-4; end
if ~isfield(ncm,'eps_reg'), ncm.eps_reg = 1e-8; end


%% RECORDERS
x_hist      = zeros(num_x, num_t);
x_non_hist  = zeros(num_x, num_t);
xd_hist     = zeros(num_x, num_t);
u_hist      = zeros(num_u, num_t);
ud_hist     = zeros(num_u, num_t);
th_ref_hist = zeros(1, num_t);

X_hist        = zeros(1, num_t);
optDone_hist  = zeros(1, num_t);
M_hist        = zeros(1, num_t);
mu_hist       = zeros(1, num_t);

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n");
tic
for t_idx = 1:num_t

    % desired input from BSC (with guard & saturation inside)
    ud = BSC_input(xd, xd_func(t(t_idx)), dxd_func(t(t_idx)), ddxd_func(t(t_idx)), param);
    % BSC 포화(추가 보강)
    ud = max(min(ud,  I_MAX), -I_MAX);

    % --- control update at controller rate ---
    if mod(t_idx-1, ctrl_stride) == 0
        e = x - xd;

        % NCM update
        try
            ncm = NCM_ctrl3(ncm, x, xd, ud, param)
        catch ME
            warning('NCM_ctrl crashed at t=%.6f: %s', t(t_idx), ME.message);
            ncm.optDone = -1;
        end

        % input matrix B
        th  = x(1);
        J   = param.J;   P = param.P;  Phi = param.Phi;
        B   = [0 0; -(3*P*Phi)/(2*J)*sin(P*th) (3*P*Phi)/(2*J)*cos(P*th)];

        % control law (no inv, with regularization & guard)
        mu_reg = max(ncm.mu, ncm.mu_min);
        W_reg  = ncm.W_bar + ncm.eps_reg*eye(size(ncm.W_bar));


            K = ncm.inv_R * B' * ((W_reg/mu_reg) \ e);
            K     = ncm.kc * K;   
       
        
            u_cmd = ud - K;


% rate limit (액츄에이터의 스펙)
du = u_cmd - u_prev;
du = max(min(du, du_max), -du_max);
u  = u_prev + du;


% hard saturation
% u  = max(min(u, I_MAX), -I_MAX);

u_prev = u;







        
        fprintf("\tControl Decision at t = %.4f, flag: %d\n", t(t_idx), ncm.optDone);
    end

    % --- record ---
    x_hist(:, t_idx)     = x;
    x_non_hist(:, t_idx) = x_non;
    xd_hist(:, t_idx)    = xd;
    u_hist(:, t_idx)     = u;
    ud_hist(:, t_idx)    = ud;

    X_hist(t_idx)        = ncm.X;
    optDone_hist(t_idx)  = ncm.optDone;
    mu_hist(t_idx)       = ncm.mu;
    th_ref_hist(t_idx)   = xd_func(t(t_idx));

    % metric condition (inv 대신 cond)
    mu_reg = max(ncm.mu, ncm.mu_min);
    W_reg  = ncm.W_bar + ncm.eps_reg*eye(size(ncm.W_bar));
    M_hist(t_idx) = cond(W_reg/mu_reg);

    % --- disturbance (smooth sign -> tanh) ---
    d_MAX  = 0.05;
    eps_s  = 1e-3;                 % smoothing scale
    d_func = @(x2) d_MAX * tanh(x2/eps_s);

    % --- integrate (Euler) ---
    x     = system_step(dt, x,     u,  d_func(x(2)),     param);
    x_non = system_step(dt, x_non, ud, d_func(x_non(2)), param);
    xd    = system_step(dt, xd,    ud, 0,                param);

    % 상태 유한성 검사 + 필요시 중단/롤백
    if any(~isfinite([x; xd; x_non]))
        warning('Non-finite state at t=%.6f. Stopping.', t(t_idx));
        if t_idx > 1
            x     = x_hist(:, t_idx-1);
            xd    = xd_hist(:, t_idx-1);
            x_non = x_non_hist(:, t_idx-1);
        end
        break;
    end

    % --- report ---
    if mod(t_idx-1, rpt_stride) == 0
        fprintf('Simulation Time: %.4f\n', t(t_idx));
    end
end
toc
fprintf("SIMULATION is Terminated\n");

%% RESULT SAVE
whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    fprintf("\nRESULT SAVING...\n");
    saveName = "results/" + whatTimeIsIt + ".mat";
    save(saveName, 'x_hist','u_hist','xd_hist','t','mu_hist','M_hist','X_hist','optDone_hist');
    fprintf("RESULT Saved: %s\n", saveName);
end

%% PLOT
if FIGURE_PLOT_FLAG
    fprintf("\nFIGURE PLOTTING...\n");
    plotter2       % 기존 플로터 사용
    fprintf("FIGURE PLOTTING Done\n");

    if FIGURE_SAVE_FLAG
        fprintf("\nFIGURE SAVING...\n");
        saveName = "figures/" + whatTimeIsIt;
        [~,~] = mkdir(saveName);
        for idx = 1:4
            f_name = saveName + "/Fig" + string(idx);
            if isvalid(figure(idx))
                saveas(figure(idx), f_name + ".png");
                exportgraphics(figure(idx), f_name + ".eps");
            end
        end
        fprintf("FIGURES Saved in %s\n", saveName);
    end
end

beep()

%% =================== LOCAL FUNCTIONS ===================

function x = system_step(dt, x, u, d, param)
    % 상태: x = [theta; omega]
    th = x(1);  w = x(2);
    ia = u(1);  ib = u(2);

    J   = param.J;   Phi = param.Phi;  P = param.P;  fv = param.fv;

    trq = -(3/2)*P*Phi*sin(P*th)*ia + (3/2)*P*Phi*cos(P*th)*ib;

    grad = [ w;
            (1/J)*(trq - fv*w - d) ];

    x = x + grad*dt;
end

function ud = BSC_input(x, r1, dr1, ddr1, param)
    % BSC (Backstepping-like) 입력계산 + 가드/포화
    k1 = 80;%1e2;  % position gain
    k2 = 40;%5e1;  % velocity gain
    % k1 = 40;  % position gain
    % k2 = 15;  % velocity gain

    x1 = x(1); x2 = x(2);
    J   = param.J;  Phi = param.Phi;  P = param.P;  fv = param.fv;

    % 유한성 가드
    if ~all(isfinite([x1,x2,r1,dr1,ddr1,J,Phi,P,fv]))
        ud = [0;0]; return
    end

    alp = dr1 - k1*(x1 - r1);
    trq = J*( -k2*(x2 - alp) - (x1 - r1) + fv*x2 + ddr1 - k1*(x2 - dr1) );

    denom = (P*Phi);
    if ~isfinite(denom) || abs(denom) < 1e-9
        ud = [0;0]; return
    end

    ud = [ -(2/3)*trq * sin(P*x1) / denom;
            +(2/3)*trq * cos(P*x1) / denom ];

    % 포화 (메인에서 I_MAX = 30 일관성)
    I_MAX = 2;
    ud = max(min(ud, I_MAX), -I_MAX);
end
