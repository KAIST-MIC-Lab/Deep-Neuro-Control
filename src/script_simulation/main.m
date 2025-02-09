% ============================================
% I understand that this is not familiar to you at the very first time.
% But, I am sure that you will get used to it soon.
% 
% ============================================

%% FASTEN YOUR SEATBELT
clear

RESULT_SAVE_FLAG = 1;
FIGURE_PLOT_FLAG = 1;
FIGURE_SAVE_FLAG = 1;

%% SIMULATION SETTING
T = 10;                 % simulation time
ctrl_dt = 1e-4;         % controller sampling time
dt = ctrl_dt * 1;       % simulation sampling time
rpt_dt = 1;             % report time (on console)
t = 0:dt:T;             % time vector

%% SYSTEM AND REFERENCE DEFINITION
x = [0;0];             % initial state
u = [0;0];             % initial input 

grad = @system_grad;

ref = @(t) [
    sin(t);
    cos(t);
];

num_x = length(x);
num_u = length(u);
num_t = length(t);

%% CONTROLLER LOAD
K = diag([2; 3]);

%% RECORDER SETTING
x_hist = zeros(num_x, num_t);
u_hist = zeros(num_u, num_t);
r_hist = zeros(num_x, num_t);

%% MAIN LOOP
for t_idx = 1:1:num_t
    r = ref(t(t_idx));
    e = x - r;

    % Control Decision
    u = -K'*e;
    
    % Record
    x_hist(:, t_idx) = x;
    u_hist(:, t_idx) = u;
    r_hist(:, t_idx) = r;

    % Step forward
    x = x + grad(x, u) * dt;

    % Report
    if mod(t_idx, rpt_dt/dt) == 0
        fprintf('Time: %.2f\n', t(t_idx))
    end
end

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG

    saveName = "results/"+whatTimeIsIt+".mat";
    save(saveName, 'x_hist', 'u_hist', 'r_hist', 't')
end

if FIGURE_PLOT_FLAG
    plotter

    if FIGURE_SAVE_FLAG
        [~,~] = mkdir("figures/"+whatTimeIsIt);

        for idx = 1:1:4   
            f_name = "figures/" + whatTimeIsIt + "/Fig" + string(idx);
    
            saveas(figure(idx), f_name + ".png")
            exportgraphics(figure(idx), f_name+'.eps')
        end
    end
end

beep()

%% LOCAL FUNCTIONS
function grad = system_grad(x, u)
    A = [0 1; -2 -3];      % system matrix
    B = eye(2);            % input matrix

    grad = A*x + B*u;
end
