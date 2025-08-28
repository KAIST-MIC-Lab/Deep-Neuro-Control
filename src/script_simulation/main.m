
%% PATH
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));
addpath(genpath([pwd filesep 'mosek/11.0/tools/platform/osxaarch64/bin']));
addpath(genpath([pwd filesep 'mosek/11.0/toolbox/r2022bom']));
% C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin
%% FASTEN YOUR SEATBELT
clear
RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder
FIGURE_PLOT_FLAG = 1;   % plot the result
FIGURE_SAVE_FLAG = 0;   % save the figure as .png and .eps

%% SIMULATION SETTING
T = 10;                 % simulation time
% ctrl_dt = 1e-3;         % controller sampling time
ctrl_dt = 1e-2;         % controller sampling time
% dt = 1/20e3;
dt = 1e-3;
rpt_dt = 1e-1;             % report time (on console)
t = 0:dt:T;             % time vector

%% SM PARAMETERS
param.A = [-2 1; -1 -0.5];
param.B = [1 0; 0 1];

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

%% SYSTEM AND REFERENCE DEFINITION
init_list = [
    [0;0]
];

sys = [];
for idx = 1:1:size(init_list, 2)
    sys{idx}.x = init_list(:, idx);
    sys{idx}.x_non = init_list(:, idx); % without controller (only feedforward)
    sys{idx}.u = transpose([0,0]);
end

xd = transpose([0,0]);
grad = @system_grad;

%% REFERENCE DEFINITION (will be tracked by BSC)
ud_func = @(t) [sin(t); cos(t)];

%% 
num_x = size(init_list, 1);      % number of states
num_u = length(ud_func(0));      % number of inputs
num_t = length(t);      % number of time steps

num_sample = size(init_list, 2); % number of initial conditions

%% CONTROLLER LOAD
for idx = 1:1:num_sample
    sys{idx}.ncm = NCM_init(ctrl_dt);
end

%% RECORDER SETTING
for idx = 1:1:num_sample
    sys{idx}.x_hist = zeros(num_x, num_t);   % state history 
    sys{idx}.x_non_hist = zeros(num_x, num_t); % state (with ud) history
    sys{idx}.u_hist = zeros(num_u, num_t);   % input history  
    sys{idx}.uSat_hist = zeros(num_u, num_t);   

    sys{idx}.X_hist = zeros(1, num_t);       % NCM: ratio (SS error)
    sys{idx}.optDone_hist = zeros(1, num_t); % NCM: optimization result flag
    sys{idx}.M_hist = zeros(1, num_t);       % NCM: norm of M (contraction metric)
    sys{idx}.mu_hist = zeros(1, num_t);      % NCM: upper bound of norm(M)
end

xd_hist = zeros(num_x, num_t); % desired state history (without disturbance)
ud_hist = zeros(num_u, num_t); % desired input history

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t

    % compute desired control input (BSC)
    % ud = BSC_input(xd, xd_func(t(t_idx)), dxd_func(t(t_idx)), ddxd_func(t(t_idx)), param);
    ud = ud_func(t(t_idx));

    % Control Decision
    % if t_idx==1 || rem(t(t_idx)/dt, ctrl_dt/dt) == 0
    if mod(t_idx, round(ctrl_dt/dt)) == 1
    
        for c_idx = 1:1:num_sample
            x = sys{c_idx}.x;
            x_non = sys{c_idx}.x_non;
            ncm = sys{c_idx}.ncm;

            % NCM controller call
            ncm = NCM_ctrl(ncm, x, xd, ud, param);  % controller call
            
            B = param.B;
            u = ud - ncm.inv_R*B' * inv(ncm.W_bar/ncm.mu) * (x-xd);

            max_u = 2;
            % uSat = max(min(u,max_u), -max_u);
            uSat = u;

            % update controller
            sys{c_idx}.ncm = ncm;
            sys{c_idx}.x = x;
            sys{c_idx}.x_non = x_non;
            sys{c_idx}.u = u;
            sys{c_idx}.uSat = uSat;

            % report on console
            fprintf("\tControl at t = %.4f, flag: %d\n", t(t_idx), ncm.optDone)
    

        end
    end

    xd_hist(:, t_idx) = xd;
    ud_hist(:, t_idx) = ud;

    % disturbance
    d_MAX = 1;  % maximum disturbance
    d_func = @(x2) -1 * d_MAX * sign(x2);  

    % step forward
    for c_idx = 1:1:num_sample
        % record to history
        sys{idx}.x_hist(:, t_idx) = x;
        sys{idx}.x_non_hist(:, t_idx) = x_non;
        sys{idx}.u_hist(:, t_idx) = u;
        sys{idx}.uSat_hist(:, t_idx) = uSat;
        sys{idx}.X_hist(t_idx) = ncm.X;  % controller gain
        sys{idx}.optDone_hist(t_idx) = ncm.optDone;  % optimization done flag
        sys{idx}.M_hist(:, t_idx) = norm(inv(ncm.W_bar/ncm.mu));
        sys{idx}.mu_hist(:, t_idx) = ncm.mu;
        
        sys{c_idx}.x = system_step(dt, sys{c_idx}.x, sys{c_idx}.uSat, d_func(sys{c_idx}.x(2)), param);
        sys{c_idx}.x_non = system_step(dt, sys{c_idx}.x_non, ud, d_func(sys{c_idx}.x_non(2)), param);
    end
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
    A = param.A;
    B = param.B;

    %%   
    grad = A*x + B*u + d;
    x = x + grad * dt;

end
