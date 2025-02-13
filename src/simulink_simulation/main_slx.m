% ********************************************
%
%   Hey! You finally reached here.
%
%   Most people typically use Simulink without MATLAB scripts.
%   However, using Simulink with a MATLAB script allows you to run simulations in the MATLAB environment and easily save the results in a .mat file for further analysis.
%   I created this example to illustrate how I use Simulink with MATLAB scripts.
%   You can see the following features:
%   - The Simulink simulation can be controlled by a MATLAB script.
%        (You can run the simulation even without opening the Simulink model.)
%   - The parameters in Simulink can be modified using a MATLAB script.
%   - The results of the Simulink simulation can be exported to the MATLAB workspace.
%
%   To try this example, follow these steps:
%   1. Run this script to start Simulink simulation (with RUN_FLAG = 1).
%   2. Check the results in the MATLAB workspace. (sim_result)
%   3. Save the results as a .mat file (with RESULT_SAVE_FLAG = 1).
%   4. Check the saved file in the results folder.
%   5. Run the plotter_individual.m script to plot the results.
%
%   Good luck!
%
%                               Myeongseok Ryu
%  	    				dding_98@gm.gist.ac.kr
%                                  09.Feb.2025
%
% ********************************************

%% FASTEN YOUR SEATBELT
clear

RUN_FLAG = 0;           % run the simulink simulation
RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder

slx_name = "main.slx";  % simulink file name

%% SIMULATION SETTING
T = 10;                 % simulation time
ctrl_dt = 1e-4;         % controller sampling time
dt = ctrl_dt * 1;       % simulation sampling time
t = 0:dt:T;             % time vector

%% REPORT SETTING
fprintf("\n")
fprintf("      *** SIMULATION INFORMATION ***\n")
fprintf("Termiation Time  : %.2f\n", T)
fprintf("Controller dt    : %.2e\n", ctrl_dt)
fprintf("Simulation dt    : %.2e\n", dt)
fprintf("\n")

%% INITIAL CONDITION
x = [0;0];              % initial state
u = [0;0];              % initial input 

%% CONTROLLER LOAD
K = diag([2; 3]);       % controller gain

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
