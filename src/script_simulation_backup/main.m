
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));

%% FASTEN YOUR SEATBELT
clear

RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder
FIGURE_PLOT_FLAG = 1;   % plot the result
FIGURE_SAVE_FLAG = 0;   % save the figure as .png and .eps

%% SIMULATION SETTING
T = 5e-2;                 % simulation time
ctrl_dt = 1e-3;         % controller sampling time
% dt = ctrl_dt * 1e-1;       % simulation sampling time
dt = 1/20e3;
rpt_dt = 1e-3;             % report time (on console)
t = 0:dt:T;             % time vector

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
% x = transpose([1,5,2,-2]);              % initial state 
x = transpose([.2,0,0,0]);              % initial state 
x_non = x;
u = transpose([0,0]);              % initial input
xd = transpose([0,0,0,0]);

grad = @system_grad;    % system gradient

% ud_func = @(t) [sin(20*t); cos(15*t)] * 10;
% ud_func = @(t) [sin(20*t); cos(15*t)] * 1;
ud_func = @(t) [1;1] * heaviside(t-0.02) + -.5;
trq_d_func = @(t) sin(20*t) * 0e-1;

num_x = length(x);      % number of states
num_u = length(u);      % number of inputs
num_t = length(t);      % number of time steps

%% CONTROLLER LOAD
NCM_init

%% RECORDER SETTING
x_hist = zeros(num_x, num_t);   % state history 
x_non_hist = zeros(num_x, num_t);   
xd_hist = zeros(num_x, num_t); % state derivative history
u_hist = zeros(num_u, num_t);   % input history  
ud_hist = zeros(num_u, num_t); % input derivative history
X_hist = zeros(1, num_t);   % 
optDone_hist = zeros(1, num_t); % optimization done history
M_hist = zeros(1, num_t);   % 

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t
    
    % Control Decision
    ud = ud_func(t(t_idx));  % reference input
    e = x - xd;  % error
    
    ncm = NCM_ctrl(ncm, x, xd, ud);  % controller call

    L = .66e-3;    % Inductance (mH)
    B = [0 0;1/L 0; 0 1/L];
    u = ud - ncm.inv_R * B' * inv(ncm.W) * e(2:end);  % control input
    % u = ud - ncm.inv_R * B' * inv(ncm.W) * e;  % control input

    % Record
    x_hist(:, t_idx) = x;
    x_non_hist(:, t_idx) = x_non;
    xd_hist(:, t_idx) = xd;
    u_hist(:, t_idx) = u;
    ud_hist(:, t_idx) = ud;
    X_hist(t_idx) = ncm.X;  % controller gain
    optDone_hist(t_idx) = ncm.optDone;  % optimization done flag
    M_hist(:, t_idx) = norm(inv(ncm.W));

    % Step forward
    trq_d = trq_d_func(t(t_idx));
    x = system_step(dt, x, u, trq_d);
    x_non = system_step(dt, x_non, ud, trq_d);
    xd = system_step(dt, xd, ud, 0);

    % Report
    if mod(t_idx, rpt_dt/dt) == 0
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
function x = system_step2(dt, x, u, d)
    %%
    th = x(1);  % theta
    thd = x(2);  % theta dot
    ia = x(3);  % Ia
    ib = x(4);  % Ib

    trq_d = d;
    
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
    grad = [
        thd;
        (1/J)*(trq - fv*thd + trq_d);
        (1/L)*(u(1) - R*ia + P*Phi*thd*sin(P*th));
        (1/L)*(u(2) - R*ib - P*Phi*thd*cos(P*th))
    ];
    
    x = x + grad * dt;

end

