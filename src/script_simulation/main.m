
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
T = 1;                 % simulation time
% T = 5;                 % simulation time
ctrl_dt = 1e-2;         % controller sampling time
dt = 1e-3;
% ctrl_dt = 1e-3;         % controller sampling time
% dt = 1e-4;
rpt_dt = 1e-1;          % report time (on console)
t = 0:dt:T;             % time vector

%% SM PARAMETERS
param.sig = 10;
param.rho = 28;
param.beta = 8/3;

% param.B = rand(3) + eye(3)*3;
param.B = eye(3);
% param.max_u = 1.5e2;
param.max_u = 75;

% disturbance
d_MAX = 50;  % maximum disturbance;
% d_func = @(t, x) [sin(t);-cos(t);sin(t)] * d_MAX;
% d_func = @(t, x) [1;1;1] * heaviside(t-1.5) * d_MAX;
d_func = @(t, x) [1;1;1] * (heaviside(t-1.5)-heaviside(t-1.5+dt))/dt * d_MAX;


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
% proposed 1 (effective space)
sys{1}.x = [15;10;20];
sys{1}.x_non = sys{1}.x;
sys{1}.ctrl_opt = 11;
sys{1}.u = [0;0;0];

% existing (large penalty)
sys{2}.x = [15;10;20];
sys{2}.x_non = sys{2}.x;
sys{2}.ctrl_opt = 32;
sys{2}.u = [0;0;0];
% % existing (small penalty)
% sys{3}.x = [15;10;20];
% sys{3}.x_non = sys{3}.x;
% sys{3}.ctrl_opt = 32;
% sys{3}.u = [0;0;0];

xd = [1;-3;5];
grad = @system_grad;

%% REFERENCE DEFINITION
% ud_func = @(t, xd) [sin(t); cos(t); -cos(t)]*0.1;
ud_func = @(t, x)  inv(param.B) * (- system_func(x, param) - 2*(x-[8.485;8.485;27]));

%% PASSIVE VARIABLES
num_x = length(xd);      % number of states
num_u = length(ud_func(0, xd));      % number of inputs
num_t = length(t);      % number of time steps

num_sample = length(sys); % number of initial conditions

%% SAMPLE SETTING
for idx = 1:1:num_sample
    sys{idx}.ncm = NCM_init(ctrl_dt, sys{idx}.ctrl_opt);  % NCM controller initialization

    sys{idx}.x_hist = zeros(num_x, num_t);                  % state history 
    sys{idx}.x_non_hist = zeros(num_x, num_t); % state (with ud) history
    sys{idx}.u_hist = zeros(num_u, num_t);   % input history  
    sys{idx}.uSat_hist = zeros(num_u, num_t);   

    sys{idx}.X_hist = zeros(1, num_t);       % NCM: ratio (SS error)
    sys{idx}.optDone_hist = zeros(1, num_t); % NCM: optimization result flag
    sys{idx}.M_hist = zeros(1, num_t);       % NCM: norm of M (contraction metric)
    sys{idx}.mu_hist = zeros(1, num_t);      % NCM: upper bound of norm(M)
    sys{idx}.contraction_flag_hist = zeros(1, num_t); % NCM: contraction condition flag
end

xd_hist = zeros(num_x, num_t); % desired state history (without disturbance)
ud_hist = zeros(num_u, num_t); % desired input history

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t

    % compute desired control input 
    ud = ud_func(t(t_idx), xd);

    % Control Decision
    if mod(t_idx, round(ctrl_dt/dt)) == 1
        for c_idx = 1:1:num_sample
            x = sys{c_idx}.x;
            x_non = sys{c_idx}.x_non;
            ncm = sys{c_idx}.ncm;

            % NCM controller call
            ncm = NCM_ctrl(ncm, x, xd, ud, param);  % controller call
            
            B = param.B;
            % u = ud - ncm.inv_R*B' * inv(ncm.W_bar/ncm.mu) * (x-xd);
            u = ud - ncm.inv_R*B' * ncm.M * (x-xd);

            max_u = param.max_u;
            % uSat = max(min(u,max_u), -max_u);
            if norm(u) > max_u; uSat = u/norm(u)*max_u; else; uSat = u; end
            % uSat = u;

            % update controller
            sys{c_idx}.ncm = ncm;
            sys{c_idx}.u = u;
            sys{c_idx}.uSat = uSat;

            % report on console
            % fprintf("\tControl at t = %.4f, flag: %d\n", t(t_idx), ncm.optDone)
        end
    end

    % step forward
    for c_idx = 1:1:num_sample
        % record to history
        sys{c_idx}.x_hist(:, t_idx) = sys{c_idx}.x;
        sys{c_idx}.x_non_hist(:, t_idx) = sys{c_idx}.x_non;
        sys{c_idx}.u_hist(:, t_idx) = sys{c_idx}.u;
        sys{c_idx}.uSat_hist(:, t_idx) = sys{c_idx}.uSat;

        sys{c_idx}.X_hist(t_idx) = sys{c_idx}.ncm.X;  % controller gain
        sys{c_idx}.optDone_hist(t_idx) = sys{c_idx}.ncm.optDone;  % optimization done flag
        sys{c_idx}.M_hist(:, t_idx) = norm(sys{c_idx}.ncm.M);
        if isfield(sys{c_idx}.ncm, 'mu')
            sys{c_idx}.mu_hist(:, t_idx) = sys{c_idx}.ncm.mu;
        end
        sys{c_idx}.contraction_flag_hist(:, t_idx) = sys{c_idx}.ncm.contraction_flag;

        % step forward
        sys{c_idx}.x = system_step(dt, sys{c_idx}.x, sys{c_idx}.uSat, d_func(t(t_idx), sys{c_idx}.x), param);
        sys{c_idx}.x_non = system_step(dt, sys{c_idx}.x_non, ud, d_func(t(t_idx), sys{c_idx}.x_non), param);
    end

    % record desired history
    xd_hist(:, t_idx) = xd;
    ud_hist(:, t_idx) = ud;

    % step forward desired state
    xd = system_step(dt, xd, ud, 0, param);

    % report on console
    if mod(t_idx, round(rpt_dt/dt)) == 1
        fprintf('Simulation Time: %.3f (%.2f %%)\n', t(t_idx), t_idx/num_t*100);
    end
end

fprintf("SIMULATION is Terminated\n")

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    fprintf("\n")
    fprintf("RESULT SAVING...\n")

    saveName = "results/"+whatTimeIsIt+".mat";
    save(saveName, 'sys', 'param', 't', 'xd_hist', 'ud_hist')

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
    % Compute system dynamics
    f = system_func(x, param);
    B = param.B;

    %%   
    grad = f + B*u + d;
    x = x + grad * dt;

end

function f = system_func(x, param)
    lr_x = x(1); lr_y = x(2); lr_z = x(3);
    sig = param.sig; rho = param.rho; beta = param.beta;

    %%
    f = [
        sig*(lr_y - lr_x);
        lr_x*(rho-lr_z)-lr_y;
        lr_x*lr_y-beta*lr_z
        ];
end