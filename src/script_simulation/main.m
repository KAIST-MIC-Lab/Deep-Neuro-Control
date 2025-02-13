% ********************************************
%
%   Hallo?,
% 
%   I understand that this may not be familiar to you at first.
%   But I’m sure you’ll get used to it soon.
% 
%   Run this script first, then explore plotter.m to see how the results are plotted.
%   Enjoy!
%
%                               Myeongseok Ryu
%  	    				dding_98@gm.gist.ac.kr
%                                  09.Feb.2025
%
% ********************************************

%% FASTEN YOUR SEATBELT
clear

RESULT_SAVE_FLAG = 1;   % save the result as a .mat file in the results folder
FIGURE_PLOT_FLAG = 1;   % plot the result
FIGURE_SAVE_FLAG = 1;   % save the figure as .png and .eps

%% SIMULATION SETTING
T = 10;                 % simulation time
ctrl_dt = 1e-4;         % controller sampling time
dt = ctrl_dt * 1;       % simulation sampling time
rpt_dt = 1;             % report time (on console)
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
fprintf("\n")

%% SYSTEM AND REFERENCE DEFINITION
x = [0;0];              % initial state
u = [0;0];              % initial input 

grad = @system_grad;    % system gradient

ref = @(t) [            % reference function
    sin(t);
    cos(t);
];

num_x = length(x);      % number of states
num_u = length(u);      % number of inputs
num_t = length(t);      % number of time steps

%% CONTROLLER LOAD
K = diag([2; 3]);       % controller gain

%% RECORDER SETTING
x_hist = zeros(num_x, num_t);   % state history 
u_hist = zeros(num_u, num_t);   % input history
r_hist = zeros(num_x, num_t);   % reference history

%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t
    % Error Calculation
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
        fprintf('Simulation Time: %.2f\n', t(t_idx))
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

%% LOCAL FUNCTIONS
function grad = system_grad(x, u)
    A = [0 1; -2 -3];      % system matrix
    B = eye(2);            % input matrix

    grad = A*x + B*u;
end
