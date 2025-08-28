%% FASTEN YOUR SEATBELT
clear


addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));

%%
RUN_FLAG = 0;           % run the simulink simulation
RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder

slx_name = "main.slx";  % simulink file name

%% SIMULATION SETTING
T = 10;                 % simulation time
ctrl_dt = 1/1e3;         % controller sampling time
dt = 1/20e3;       % simulation sampling time
t = 0:dt:T;             % time vector

%% REPORT SETTING
fprintf("\n")
fprintf("      *** SIMULATION INFORMATION ***\n")
fprintf("Termiation Time  : %.2f\n", T)
fprintf("Controller dt    : %.2e\n", ctrl_dt)
fprintf("Simulation dt    : %.2e\n", dt)
fprintf("\n")

%% INITIAL CONDITION
x = [0;0;0;0];              % initial state
u = [0;0];              % initial input 

%% CONTROLLER
x_num = 3;
W_init = zeros(x_num);
X_init = 1;

%% MAIN SIMULATION RUN
if RUN_FLAG
    fprintf("SIMULINK SIMULATION is Running...\n")

    sim_result = sim(slx_name);
    
    fprintf("SIMULINK SIMULATION is Done!\n")

    beep()
end

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    fprintf("\n")
    fprintf("RESULT SAVING...\n")

    saveName = "results/"+whatTimeIsIt+".mat";

    logsout = sim_result.logsout;
    save(saveName, "logsout", "T");

    fprintf("RESULT SAVED as %s\n", saveName)
end

beep()
