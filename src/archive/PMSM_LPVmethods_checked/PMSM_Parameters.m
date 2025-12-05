%% PMSM General model

% Initialization
theta0 = 0;

% Parameters
Nr = 4000;         % pole numbers
J = 2.29e-3;    % Inertia
B = 0.0001  ;      % Viscous friction coefficient
Km = 0.18;     % Motor torque constant
R = 0.8;        % Phase resistance
L = 1.16e-1;    % Phase inductance

% Control parameters
Vmax = 50;      % Openloop/Conventional microstepping
ControlMode_Mech = 1;       % Mechanical control mode: [0] TorqueControl // [1] CurrentControl
ControlMode_Elec = 2;       % Electrical control mode: [0] OpenLoop // [1] PID current control with microstepping voltage // [2] Algorithm
Current_Control = 0;        % [0] FOC // [1] FWC
OmegaStar_Selector = 0; 

%% Controller gains

% Current PID Control
kp_elec = 10;   % Current PID control P gain
ki_elec = 0.01; % Current PID control I gain
kd_elec = 0.05; % Current PID control D gain

% CST13'
k1_CST = 0;     % omega star k1 gain
k2_CST = 0;     % omega star k2 gain
k3_CST = 0;     % omega star k3 gain

