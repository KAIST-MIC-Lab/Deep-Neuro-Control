clear

%% SIMULATION SETTING
T = 1;                 % simulation time
ctrl_dt = 1e-4;         % controller sampling time
dt = 1e-5;
rpt_dt = 1e-1;          % report time (on console)
t = 0:dt:T;             % time vector

%%
x = [15;10;20];
u0 = [0;0;0];
u = u0;

%%
xd_func = @(t) [sin(t); cos(t); tanh(t)] .* [20;16;12];
xd_dot_func = @(t) [cos(t); -sin(t); sech(t)^2] .* [20;16;12];

%%
param.sig = 10;
param.rho = 28;
param.beta = 8/3;

u_max = 160;
sat_func = @(u, u_max) max(min(u, u_max), -u_max);

function d = d_func(dt, t, x)
    if t>=0.5 && t<0.5+dt
        d = [1;1;1]*10/dt;
    else
        d = [0;0;0];
    end
end

%%
epsilon = 1e-3;
K = diag([1, 1, 1]) * 1e1;

%%
num_t = length(t);
num_x = size(x, 1);
num_u = size(u, 1);

%% 
x_hist = zeros(num_x, num_t); x_hist(:,1) = x;
u_hist = zeros(num_u, num_t); u_hist(:,1) = u;
xd_hist = zeros(num_x, num_t); xd_hist(:,1) = xd_func(t(1));

%% 
for t_idx = 1:1:num_t

    xd = xd_func(t(t_idx));
    [f, B] = system_dynamics(x, u, param);

    % control update
    if mod(t_idx, round(ctrl_dt/dt)) == 1
        xd_dot = xd_dot_func(t(t_idx));

        grad_u = -B' * (K*(x-xd) + f - xd_dot + B*sat_func(u, u_max));
        u = u + (1/epsilon*grad_u) * ctrl_dt;
    end

    % system update
    d = d_func(dt, t(t_idx), x);
    x = x + (f+B*sat_func(u, u_max) + d) * dt;

    % record
    x_hist(:, t_idx) = x;
    u_hist(:, t_idx) = u;
    xd_hist(:, t_idx) = xd;
end

%%
figure(1); clf; hold on; box on; grid on;
tl = tiledlayout(3,1);
tl.Padding = 'compact';
tl.TileSpacing = 'compact';

nexttile;
plot(t, x_hist(1,:), 'LineWidth', 1.5); hold on
plot(t, xd_hist(1,:), '--', 'LineWidth', 1.5);
ylabel('x_1');
legend('x_1', 'x_{d1}');
title('State Tracking');
set(gca, 'FontSize', 12);

nexttile;
plot(t, x_hist(2,:), 'LineWidth', 1.5); hold on
plot(t, xd_hist(2,:), '--', 'LineWidth', 1.5);
ylabel('x_2');
legend('x_2', 'x_{d2}');
set(gca, 'FontSize', 12);

nexttile;
plot(t, x_hist(3,:), 'LineWidth', 1.5); hold on
plot(t, xd_hist(3,:), '--', 'LineWidth', 1.5);
ylabel('x_3');
xlabel('Time (s)');
legend('x_3', 'x_{d3}');
set(gca, 'FontSize', 12);   

figure(2); clf; hold on; box on; grid on;
plot(t, u_hist(1,:), 'LineWidth', 1.5); hold on
plot(t, u_hist(2,:), 'LineWidth', 1.5);
plot(t, u_hist(3,:), 'LineWidth', 1.5);
yline(u_max, '--k', 'u_{max}');
yline(-u_max, '--k', '-u_{max}');
ylabel('Control Input');
xlabel('Time (s)');
legend('u_1', 'u_2', 'u_3', 'Location', 'best');
set(gca, 'FontSize', 12);
ylim([-u_max u_max]*1.2);

%% 
function [f, B] = system_dynamics(x, u, param)
    f = system_func(x, param);
    B = eye(3) * 1e0;
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