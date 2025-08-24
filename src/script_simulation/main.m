
%% PATH
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));
addpath(genpath([pwd filesep 'mosek/11.0/tools/platform/osxaarch64/bin']));
addpath(genpath([pwd filesep 'mosek/11.0/toolbox/r2022bom']));

%% FASTEN YOUR SEATBELT
clear
RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder
FIGURE_PLOT_FLAG = 1;   % plot the result
FIGURE_SAVE_FLAG = 0;   % save the figure as .png and .eps

%% SIMULATION SETTING
T = 3;                 % simulation time
% ctrl_dt = 1e-3;         % controller sampling time
ctrl_dt = 1e-2;         % controller sampling time
% dt = 1/20e3;
dt = 1e-6;
rpt_dt = 1e-1;             % report time (on console)
t = 0:dt:T;             % time vector

%% SM PARAMETERS
param.L = .66e-3;    % Inductance (mH)
param.R = 0.251;     % Resistance (Ohm)
param.J = 3.24e-5;   % Inertia (kg.m^2)
param.Phi = 16.8e-3; % Flux (Wb)
param.P = 4;         % Pole pairs
param.fv = 2e-3;     % Viscous friction (N.m.s/rad)
    

%% REPORT SETTING
fprintf("\n")
fprintf("      *** SIMULATION INFORMATION ***\n")
fprintf("Termiation Time  : %.2f\n", T)
fprintf("Controller dt    : %.2e\n", ctrl_dt)
fprintf("Simulation dt    : %.2e\n", dt)
fprintf("Report dt        : %.2e\n", rpt_dt)
fprintf("\n")
fprintf("RESULT_SAVE_FLAG : %d\n", RESULT_SAVE_FLAG)
fprintf("FIGURE_PLOT_FLAG : %d\n", FIGURE_PLOT_FLAG)
fprintf("FIGURE_SAVE_FLAG : %d\n", FIGURE_SAVE_FLAG)
% fprintf("\n")

%% SYSTEM AND REFERENCE DEFINITION
x = transpose([0,0]);              % initial state 
x_non = x;
u = transpose([0,0]);              % initial input
xd = transpose([0,0]);

grad = @system_grad;    % system gradient

%% REFERENCE DEFINITION (will be tracked by BSC)
xd_func = @(t) 2*(cos(t)-1) + 1*(cos(10*t)-1);
dxd_func = @(t) -2*sin(t) - 1*10*sin(10*t);
ddxd_func = @(t) -2*cos(t) - 1*100*cos(10*t);

%% 
num_x = length(x);      % number of states
num_u = length(u);      % number of inputs
num_t = length(t);      % number of time steps

%% CONTROLLER LOAD
NCM_init

%% RECORDER SETTING
x_hist = zeros(num_x, num_t);   % state history 
x_non_hist = zeros(num_x, num_t); % state (with ud) history
xd_hist = zeros(num_x, num_t); % desired state history (without disturbance)
u_hist = zeros(num_u, num_t);   % input history  
ud_hist = zeros(num_u, num_t); % desried input history
th_ref_hist = zeros(1, num_t); % reference angle history (will be tracked by BSC)

X_hist = zeros(1, num_t);       % NCM: raxio (SS error)
optDone_hist = zeros(1, num_t); % NCM: optimization result flag
M_hist = zeros(1, num_t);       % NCM: norm of M (contraction metric)
mu_hist = zeros(1, num_t);      % NCM: upper bound of norm(M)

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t

    % compute desired control input (BSC)
    ud = BSC_input(xd, xd_func(t(t_idx)), dxd_func(t(t_idx)), ddxd_func(t(t_idx)), param);

    % Control Decision
    % if t_idx==1 || rem(t(t_idx)/dt, ctrl_dt/dt) == 0
    if mod(t_idx, round(ctrl_dt/dt)) == 1
    
        % compute error
        e = x - xd;  % error
        
        % NCM controller call
        ncm = NCM_ctrl(ncm, x, xd, ud, param);  % controller call
        
        % compute input matrix B
        th = x(1);
        J = param.J;    % Inductance (mH)
        P = param.P;         % Pole pairs
        Phi = param.Phi; % Flux (Wb)
        B = [0 0; -(3*P*Phi)/(2*J)*sin(P*th) (3*P*Phi)/(2*J)*cos(P*th)];  % input matrix
        
        % compute control input
        u = ud - ncm.inv_R*B' * inv(ncm.W_bar/ncm.mu) * e; 
        % u = ud - ncm.inv_R*B' * ((ncm.W_bar/ncm.mu) \ e);

        % report on console
        fprintf("\tControl Decision at t = %.4f, flag: %d\n", t(t_idx), ncm.optDone)

    end

    % record to history
    x_hist(:, t_idx) = x;
    x_non_hist(:, t_idx) = x_non;
    xd_hist(:, t_idx) = xd;
    u_hist(:, t_idx) = u;
    ud_hist(:, t_idx) = ud;
    X_hist(t_idx) = ncm.X;  % controller gain
    optDone_hist(t_idx) = ncm.optDone;  % optimization done flag
    M_hist(:, t_idx) = norm(inv(ncm.W_bar/ncm.mu));
    mu_hist(:, t_idx) = ncm.mu;
    th_ref_hist(t_idx) = xd_func(t(t_idx));  % reference angle

    % disturbance
    d_MAX = 0.01;  % maximum disturbance
    d_func = @(x2) d_MAX * sign(x2);  

    % step forward
    x = system_step(dt, x, u, d_func(x(2)), param);
    x_non = system_step(dt, x_non, ud, d_func(x_non(2)), param);
    xd = system_step(dt, xd, ud, 0, param);
    

    % report on console
    if mod(t_idx, round(rpt_dt/dt)) == 1
        fprintf('Simulation Time: %.4f\n', t(t_idx))
    end
end

fprintf("SIMULATION is Terminated\n")

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    fprintf("\n")
    fprintf("RESULT SAVING...\n")

    saveName = "results/"+whatTimeIsIt+".mat";
    save(saveName, 'x_hist', 'u_hist', 'r_hist', 't')

    fprintf("RESULT is Saved as \n \t%s\n", saveName)
end

if FIGURE_PLOT_FLAG
    fprintf("\n")
    fprintf("FIGURE PLOTTING...\n")

    plotter

    fprintf("FIGURE PLOTTING is Done\n")

    if FIGURE_SAVE_FLAG
        fprintf("\n")
        fprintf("FIGURE SAVING...\n")
        
        saveName = "figures/"+whatTimeIsIt;
        [~,~] = mkdir(saveName);

        for idx = 1:1:4   
            f_name = saveName + "/Fig" + string(idx);
    
            saveas(figure(idx), f_name + ".png")
            exportgraphics(figure(idx), f_name+'.eps')
        end

        fprintf("FIGURE is Saved in \n \t%s\n", saveName)
    end
end

beep()

%% LOCAL FUNCTION
function x = system_step(dt, x, u, d, param)
    %%
    th = x(1);  % theta
    w = x(2);  % theta dot

    ia = u(1);  % Ia
    ib = u(2);  % Ib

    %% 
    trq_d = d;

    %%
    % L = param.L;    % Inductance (mH)
    % R = param.R;     % Resistance (Ohm)
    J = param.J;   % Inertia (kg.m^2)
    Phi = param.Phi; % Flux (Wb)
    P = param.P;         % Pole pairs
    fv = param.fv;     % Viscous friction (N.m.s/rad)
    
    %% 
    trq = -(3/2)*P*Phi*sin(P*th)*ia + (3/2)*P*Phi*cos(P*th)*ib;

    %%
    grad = [
        w;
        (1/J)*(trq - fv*w - trq_d);
    ];
    
    x = x + grad * dt;

end

function ud = BSC_input(x, xd, dxd, ddxd, param)
    %%
    k1 = 1e2;  % position gain
    k2 = 5e1;  % velocity gain
    
    %%
    x1 = x(1);  % theta
    x2 = x(2);  % theta dot

    r1 = xd(1);  

    dr1 = dxd(1);  
    ddr1 = ddxd(1);  

    %%
    % L = param.L;    % Inductance (mH)
    % R = param.R;     % Resistance (Ohm)
    J = param.J;   % Inertia (kg.m^2)
    Phi = param.Phi; % Flux (Wb)
    P = param.P;         % Pole pairs
    fv = param.fv;     % Viscous friction (N.m.s/rad)

    %%
    alp = dr1 - k1 * (x1-r1);

    trq = J * ( ...
        -k2 * (x2 - alp) ... % - (x1 - r1) ...
        + fv*x2 ... %+ ddr1 - k1*(x2 - dr1) ...
    );

    ud = [
        -(2/3)*trq * sin(P*x1) / (P*Phi);
        +(2/3)*trq * cos(P*x1) / (P*Phi)
    ];
end