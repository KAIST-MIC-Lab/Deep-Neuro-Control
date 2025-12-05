
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

PAPER_FIGURE_PLOT_FLAG = 0;   % plot the result

%% SIMULATION SETTING
T = .5;                 % simulation time
ctrl_dt = 1/200;         % controller sampling time
dt = ctrl_dt/1e1;
rpt_dt = 1e-1;          % report time )(on console)
t = 0:dt:T;             % time vector

%% SM PARAMETERS
param.sig = 10;
param.rho = 28;
param.beta = 8/3;

param.B = eye(3);
param.max_u = 200;

% disturbance
d_MAX = 1000;  % maximum disturbance;
% d_func = @(t, x) [1;1;1] * (heaviside(t,.1)-heaviside(t,.1+dt)) * d_MAX;
d_func = @(t, x) [1;1;1] * (heaviside(t,.2)-heaviside(t,.2+20*dt)) * d_MAX;
% d_func = @(t, x) randn(3,1) * heaviside(t,.1) * d_MAX;

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

%% INITIAL POINTS SETTING
% x = [10;5;13];     % initial state
x = [10;10;10];     % initial state
xd = [11;4;15];             % desired initial state
u = [0;0;0];        % initial input

%% PASSIVE CONSTANTS
num_x = length(x); % number of states
num_u = length(u); % number of inputs
num_t = length(t); % number of time steps

%% TEST CASE SETTING
% ctrl_opt 11: proposed 1 (effective space)
% ctrl_opt 32: existing (small penalty)
test_cases = [
    % struct('name', 'c1', 'ctrl_opt', 31, 'dist_on', 0, 'color', [.5 .5 .5]), ...
    struct('name', 'C3', 'ctrl_opt', 0, 'dist_on', 1, 'color', [.5 .5 .5]), ...
    struct('name', 'C2', 'ctrl_opt', 32, 'dist_on', 1, 'color', [0 1 1]), ...
    struct('name', 'C1', 'ctrl_opt', 11, 'dist_on', 1, 'color', [0 0 1]), ...
    ];
case_num = length(test_cases);

for c_idx = 1:case_num
    recs{c_idx} = struct();
    for field_idx = 1:length(fieldnames(test_cases(c_idx)))
        field_name = string(fieldnames(test_cases(c_idx)));
        recs{c_idx}.(field_name(field_idx)) = test_cases(c_idx).(field_name(field_idx));
    end

    %% STATE
    recs{c_idx}.x = x;
    recs{c_idx}.x_non = recs{c_idx}.x;
    recs{c_idx}.u = u;

    recs{c_idx}.ncm = NCM_init(ctrl_dt, recs{c_idx}.ctrl_opt);

    %% HISTORY RECORDING SETTING
    recs{c_idx}.x_hist = zeros(num_x, num_t); % state history
    recs{c_idx}.x_non_hist = zeros(num_x, num_t); % state (with ud) history
    recs{c_idx}.u_hist = zeros(num_u, num_t); % input history
    recs{c_idx}.uSat_hist = zeros(num_u, num_t);

    recs{c_idx}.x_hist(:,1) = x;
    recs{c_idx}.x_non_hist(:,1) = x;

    recs{c_idx}.X_hist = zeros(1, num_t); % NCM: ratio (SS error)
    recs{c_idx}.optDone_hist = zeros(1, num_t); % NCM: optimization result flag
    recs{c_idx}.M_hist = zeros(1, num_t); % NCM: norm of M (contraction metric)
    recs{c_idx}.nu_hist = zeros(1, num_t); % NCM: upper bound of norm(M)
    recs{c_idx}.contraction_flag_hist = zeros(1, num_t); % NCM: contraction condition flag

    recs{c_idx}.cmp_t_hist = zeros(1, num_t); % computation time history

end

%% DESIRED TRAJECTORY AND CONTROL INPUT FUNCTION
% ud_func = @(t, xd) [sin(t); cos(t); -cos(t)]*0.1;
ud_func = @(t, x)  param.B \ (- system_func(x, param) - 2*(x-[8.485;8.485;27]));

xd_hist = zeros(num_x, num_t); % desired state history (without disturbance)
ud_hist = zeros(num_u, num_t); % desired input history

%% MAIN LOOP
fprintf("\n")
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t

    % control decision
    if mod(t_idx, round(ctrl_dt/dt)) == 1
        ud = ud_func(t(t_idx), xd);
        
        for c_idx = 1:1:case_num
            x = recs{c_idx}.x; x_non = recs{c_idx}.x_non;
            ncm = recs{c_idx}.ncm;

            % NCM controller call
            ncm = NCM_ctrl(ncm, x, xd, ud, param);  % controller call
            
            B = param.B;
            % u = ud - ncm.inv_R*B' * inv(ncm.W_bar/ncm.nu) * (x-xd);
            % u = ud - ncm.inv_R*B' * ncm.M * (x-xd);
            u = ncm.u;

            max_u = param.max_u;
            % uSat = max(min(u,max_u), -max_u);
            if norm(u) > max_u
                uSat = u/norm(u)*max_u; 
            else 
                uSat = u;
            end
            % uSat = u;

            % update controller
            recs{c_idx}.ncm = ncm;
            recs{c_idx}.u = u;
            recs{c_idx}.uSat = uSat;

            % report on console
            % fprintf("\tControl at t = %.4f, flag: %d\n", t(t_idx), ncm.optDone)
        end
    end

    % step forward
    for c_idx = 1:1:case_num
        % record to history
        recs{c_idx}.x_hist(:, t_idx) = recs{c_idx}.x;
        % recs{c_idx}.x_non_hist(:, t_idx) = recs{c_idx}.x_non;
        recs{c_idx}.u_hist(:, t_idx) = recs{c_idx}.u;
        recs{c_idx}.uSat_hist(:, t_idx) = recs{c_idx}.uSat;

        recs{c_idx}.X_hist(t_idx) = recs{c_idx}.ncm.X;  % controller gain
        recs{c_idx}.optDone_hist(t_idx) = recs{c_idx}.ncm.optDone;  % optimization done flag
        recs{c_idx}.M_hist(:, t_idx) = norm(recs{c_idx}.ncm.M);
        if isfield(recs{c_idx}.ncm, 'nu')
            recs{c_idx}.nu_hist(:, t_idx) = recs{c_idx}.ncm.nu;
        end
        recs{c_idx}.contraction_flag_hist(:, t_idx) = recs{c_idx}.ncm.contraction_flag;
        
        recs{c_idx}.cmp_t_hist(:, t_idx) = recs{c_idx}.ncm.cmp_t;  % computation time

        % step forward
        if recs{c_idx}.dist_on; d = d_func(t(t_idx), recs{c_idx}.x); else; d = zeros(x_num); end
        recs{c_idx}.x = system_step(dt, recs{c_idx}.x, recs{c_idx}.uSat, d, param);
        % recs{c_idx}.x_non = system_step(dt, recs{c_idx}.x_non, ud, d, param);

        % fprintf("t: %.3f, Case %d, x1: %.3f, x2: %.3f, x3: %.3f\n", t(t_idx)*1e3, c_idx, recs{c_idx}.x(1), recs{c_idx}.x(2), recs{c_idx}.x(3));
    end

    % step forward desired state
    xd = system_step(dt, xd, ud, 0, param);
     
    % record desired history
    xd_hist(:, t_idx) = xd;
    ud_hist(:, t_idx) = ud;

    % report on console
    if mod(t_idx, round(rpt_dt/dt)) == 1
        fprintf('Time: %.3f (%.2f %%)\n', t(t_idx), t_idx/num_t*100);
    end
end

fprintf("SIMULATION is Terminated\n")

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    fprintf("\n")
    fprintf("RESULT SAVING...\n")

    saveName = "results/"+whatTimeIsIt+".mat";
    save(saveName, 'recs', 'param', 't', 'xd_hist', 'ud_hist')

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


%% PLOT FOR PAPER
if PAPER_FIGURE_PLOT_FLAG
    plotter4paper
end

%% 
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

%% LOCAL FUNCTIONS
% Heaviside function
function y = heaviside(t, d)
    if t < d
        y = 0;
    else
        y = 1;
    end
end