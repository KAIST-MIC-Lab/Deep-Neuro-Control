%% FASTEN YOUR SEATBELT
clear

RUN_FLAG = 1;
RESULT_SAVE_FLAG = 1;

slx_name = "main.slx";

%% SIMULATION SETTING
T = 10;                 % simulation time
ctrl_dt = 1e-4;         % controller sampling time
dt = ctrl_dt * 1;       % simulation sampling time
t = 0:dt:T;             % time vector

%% INITIAL CONDITION
x = [0;0];             % initial state
u = [0;0];             % initial input 

%% CONTROLLER LOAD
K = diag([2; 3]);

%% MAIN SIMULATION RUN
if RUN_FLAG
    fprintf("SIMULINK SIMULATION RUNNING...\n")
    sim_result = sim(slx_name);
    fprintf("SIMULINK SIMULATION is Terminated\n")
    beep
end

%% RESULT REPORT AND SAVE
whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));

if RESULT_SAVE_FLAG
    saveName = "results/"+whatTimeIsIt+".mat";

    logsout = sim_result.logsout;
    save(saveName, "logsout", "T");
end

beep()
