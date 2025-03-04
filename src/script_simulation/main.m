%% FASTEN YOUR SEATBELT
clear

RESULT_SAVE_FLAG = 0;   % save the result as a .mat file in the results folder
FIGURE_PLOT_FLAG = 1;   % plot the result
FIGURE_SAVE_FLAG = 0;   % save the figure as .png and .eps

%% SIMULATION SETTING
T = 1;                 % simulation time
ctrl_dt = 1e-6;         % controller sampling time
dt = ctrl_dt * 1e0;       % simulation sampling time
rpt_dt = .1;             % report time (on console)
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
x = [0;2.5];              % initial state 
u = 0;              % initial input
Delta = [0;0];          % initial fin deflection 
mu_z = 0;               % initial normal force coefficient

grad = @system_grad;    % system gradient

ref = @(t) [            % reference function
    sin(t)+2.5;
];

num_x = length(x);      % number of states
num_u = length(u);      % number of inputs
num_t = length(t);      % number of time steps

%% CONTROLLER LOAD
K = 1e3;       % controller gain

%% RECORDER SETTING
x_hist = zeros(num_x, num_t);   % state history 
                                % [pitch rate; angle of attack]
u_hist = zeros(num_u, num_t);   % input history
r_hist = zeros(num_x, num_t);   % reference history
Del_hist = zeros(num_x, num_t); % fin deflection history
mu_hist = zeros(1, num_t);      % normal force coefficient history
%% MAIN LOOP
fprintf("SIMULATION RUNNING...\n")

for t_idx = 1:1:num_t
    % Error Calculation
    r = ref(t(t_idx));
    e = x(2) - r;

    % Control Decision
    u = -K'*e;
    
    % Record
    x_hist(:, t_idx) = x;
    u_hist(:, t_idx) = u;
    r_hist(:, t_idx) = r;
    Del_hist(:, t_idx) = Delta;
    mu_hist(t_idx) = mu_z;

    % Step forward
    [x,Delta,mu_z] = system_step(dt, x, u, Delta);

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
function [x,Delta,mu_z] = system_step(dt, x, u, Delta)
    q = x(1);                  % pitch rate [deg/s]      
    alp = x(2);                % angle of attack [deg]

    % assert(alp < 20 && alp > -20, "Angle of Attack is out of range")

    %% FIN DYNAMICS
    delta_c = u;            % commanded fin deflection [deg]
    omega_a = 150;          % actuator bandwidth [rad/s]

    Delta_grad = ([0 1; 0 -1.4*omega_a]*Delta + omega_a^2*[-1; delta_c]);
    Delta = Delta + Delta_grad * dt;
    delta = Delta(1);       % fin deflection [deg]

    %% SYSTEM PARAMETERS
    d = .75;                % reference diameter [ft]
    f = 180/pi;             % radians-to-degrees conversion
    g = 32.2;               % acceleration due to gravity [ft/s^2]
    Iyy = 182.5;            % pitch moment of inertia [slug-ft^2]
    Q = 6132.8;             % dynamic pressure [lb/ft^2]
    S = .44;                % reference area [ft^2]
    V = 3109.3;             % velocity [ft/s]
    W = 450;                % weight [lb]

    %% FORCE CALCULATE
    phi_z = 103e-6*alp^3 - 945e-5*alp*abs(alp) -.170*alp;
    phi_m = 215e-6*alp^3 -195e-4*alp*abs(alp) +.051*alp;

    bm = -.206; bz = -.034;

    Cz = phi_z + bz*delta;
    % Cm = phi_m + bm*delta;

    Z = Cz*Q*S;            % normal force [lb]
    % m = Cm*Q*S*d;           % pitch moment [ft-lb]

    mu_z = Z/W;             % normal force coefficient

    %%
    grad_tmp = [
        f*g*Q*S*cos(alp/f)/W/V;
        f*Q*S*d/Iyy
    ];
    grad = grad_tmp .* [phi_z; phi_m] ...
        + [0 1;0 0] * [alp; q] + grad_tmp .* [bz; bm] * delta;
    
    x = x + grad * dt;

end
