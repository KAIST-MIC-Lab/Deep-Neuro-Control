clear

slx_name = "main.slx";
T = .96;
seed = 1; rng(seed);

RUN_FLAG = 1;
SAVE_FLAG = 1;

%% SYSTEM PARAMETERS LOAD
% data_path = 'machine_pmsm_mtpx.mat';
data_path = 'machine_IPM01.mat';
param = load(data_path);

% current vectors / grids
isd = param.machine.psi.s.arg.isd;
isq = param.machine.psi.s.arg.isq;

% flux linkages
PSISD = param.machine.psi.s.d(:,:,1);
PSISQ = param.machine.psi.s.q(:,:,1);

% differential inductances
LSDD = param.machine.L.s.dd(:,:,1);
LSDQ = param.machine.L.s.dq(:,:,1);
LSQD = param.machine.L.s.qd(:,:,1);
LSQQ = param.machine.L.s.qq(:,:,1);

Rs = param.machine.Rs;
np = param.machine.nP;
kappa = param.machine.kappa;

%% INITIAL CONDITION
init_x = [0; 0];

%% CONTROLLER SETTING
% ctrl_dt = 1e-4;%1/8e3;
ctrl_dt = 1/8e3;

init_range = 1e-3;
NN_size = [
        % layer info (node numbers)
        4 % input layer
        8
        2 % output layer
]; 

l_size = length(NN_size); % layer number
t_size = sum(NN_size(1:end-1)); % total tape number
v_size_list = zeros(l_size-1 ,1); % layer's weight number   
for idx = 1:1:l_size-1
    v_size_list(idx) = (NN_size(idx)+1) * NN_size(idx+1);      
end
v_size = sum(v_size_list); % total weight number

lbd_size = l_size-1 + 1; % multiplier number

u_max = 30;
V_max = [4;20] * 1e0;

alp = 3e1; % learning rate
beta = [ ... % update rate of Lagrange multipliers
    1e3 ... % norm constraint of first layer's weight
    1e3 ... % norm constraint of second layer's weight
    0e+1 ... % norm constraint of control input
] * 1e0;

%% NEURAL NETWORK INITIALIZATION
th = (rand(v_size,1)-1/2)*2*init_range;
lbd = zeros(lbd_size, 1);

%% SIMULINK RUN
if RUN_FLAG
    fprintf("SIMULINK SIMULATION RUNNING...\n")
    sim_result = sim(slx_name);
    fprintf("SIMULINK SIMULATION is Terminated\n")
    beep
    
    if SAVE_FLAG
        whatTimeIsIt = string(datetime('now','Format','d-MMM-y_HH-mm-ss'));
        saveName = "sim_result/"+whatTimeIsIt+".mat";
    
        logsout = sim_result.logsout;
        save(saveName, "logsout", "u_max", "V_max", "T");
    end
end


