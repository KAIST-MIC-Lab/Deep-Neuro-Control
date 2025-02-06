clear
addpath('utils_NN')
addpath('machine_data/')

%% SIMULATION SETTING
T = 1;
% ctrl_dt = 1e-6;
% dt = ctrl_dt;
ctrl_dt = 1/8e3;
dt = ctrl_dt * 1/100;
% ctrl_dt = dt;
rpt_dt = 0.1;
t = 0:dt:T;
% seed = 1; rng(seed);

%% SYSTEM DYNAMICS
w = 0;
i = [0; 0];

% grad = @sysGrad_1st;
% grad = @sysGrad_linear;
grad = @sysGrad_lookup;
data_path = 'machine_data/machine_pmsm_mtpx.mat';
param = load(data_path);

ref = @(t) 2e2*min(1, (max(0, (t-0.1)*5)));

% ref = @(t) 1e1*(heaviside(t-1)-heaviside(t-3));
% ref = @(t) 8e2*(heaviside(t-1)-heaviside(t-3));
% ref = @(t) 1e3*(heaviside(t-1)-heaviside(t-3));

% ref = @(t) sin(10*t)*1e2;

% ref = @(t) 1e3*(heaviside(t-1)-heaviside(t-3)) ...
%     + 5e2*(heaviside(t-5)-heaviside(t-7)) ...
%     + 7e2*(heaviside(t-9)-heaviside(t-11));

%%
nnOpt1 = nn_opt_load1(ctrl_dt);
nn1 = nn_init(nnOpt1);

nnOpt2 = nn_opt_load2(ctrl_dt);
nn2 = nn_init(nnOpt2);

%%
w_hist = zeros(1, length(t));
i_hist = zeros(2, length(t));
wr_hist = zeros(1, length(t));
ir_hist = zeros(2, length(t));
v_hist = zeros(2, length(t));
vSat_hist = zeros(2, length(t));
trq_hist = zeros(1, length(t));

L1_hist = zeros(length(nnOpt1.beta), length(t));
th1_hist = zeros(nnOpt1.l_size-1, length(t));
L2_hist = zeros(length(nnOpt2.beta), length(t));
th2_hist = zeros(nnOpt2.l_size-1, length(t));
%%
LPF = 0.01; pre_wr = ref(0); pre_id = [0;0];
w_const = 0;

for t_idx = 2:1:length(t)
    r = LPF*pre_wr + (1 - LPF)*ref(t(t_idx));   

    if t_idx==2 || rem(t(t_idx)/dt, ctrl_dt/dt) == 0
        we = w-r;
        % fprintf("C on t: %.5f\n", t(t_idx));
        % x_in1 = [w;i; r];
        x_in1 = [r];
        [nn1, u_NN1, ~] = nn_forward(nn1, nnOpt1, x_in1);
        [nn1, nnOpt1, dot_L1, ~] = nn_backward(nn1, nnOpt1, w, we, u_NN1, trqCalc(i));

        id = u_NN1;

        % id = [
        %     sin(t(t_idx) * 30); 
        %     -sin(t(t_idx) * 30)
        % ] * 15e0/sqrt(2);
        
        % id = [
        %     +(+heaviside(t(t_idx)-.5)-heaviside(t(t_idx)-1)-heaviside(t(t_idx)-1.5)+heaviside(t(t_idx)-2));
        %     -(+heaviside(t(t_idx)-.5)-heaviside(t(t_idx)-1)-heaviside(t(t_idx)-1.5)+heaviside(t(t_idx)-2));
        % ] * 30;

        Del_t = 0.1;
        Del_xd = 2.5;

        k1 = fix( (t(t_idx)+Del_t/2)/Del_t );
        k2 = fix( t(t_idx)/Del_t );

        id = [
            (-1)^k1 * k1*Del_xd;
            (-1)^k2 * k2*Del_xd;
        ];


        id = LPF*pre_id + (1 - LPF)*id;   

        % id = [
        %      +min(1, max(0, (t(t_idx)-.1)*15))-min(1, max(0, (t(t_idx)-.3)*15))
        %      -min(1, max(0, (t(t_idx)-.1)*15))+min(1, max(0, (t(t_idx)-.3)*15))
        % ] * 30/sqrt(2); 


        % i = id; 
        % id = -[1;1] * we * 0.1;
        % v = [0;0];

        w = w_const;

        ie = i - id;
        x_in2 = [w;i;id];  
        [nn2, u_NN2, ~] = nn_forward(nn2, nnOpt2, x_in2);
        [nn2, nnOpt2, dot_L2, ~] = nn_backward(nn2, nnOpt2, i, ie, u_NN2, trqCalc(i));

        v = u_NN2;

        % v = -ie * 10;z

        pre_wr = r;
        pre_id = id;
    end

    % saturation
    if norm(v) > nnOpt2.cstr.u_ball
        v_sat = v/norm(v) * nnOpt2.cstr.u_ball;
    else 
        v_sat = v;
    end

    % grad_val = grad(w, i);
    % w = w + grad_val * dt;
    % % % 
   
    % grad_val = grad(w, i, v_sat);
    % w = w + grad_val(1) * dt;
    % i = i + grad_val(2:3) * dt;
    % i = id; 
    grad_val = grad(w, i, v_sat, param);
    w = w + grad_val(1) * dt;
    i = i + grad_val(2:3) * dt;
    % i = id;
    % max_err = 1e2;
    % if (norm(w) > max_err) || (norm(i) > max_err) || (norm(v) > max_err)
    %     break
    % end

    w_hist(t_idx) = w;
    i_hist(:, t_idx) = i;

    wr_hist(t_idx) = r;
    ir_hist(:, t_idx) = id;

    v_hist(:, t_idx) = v;
    vSat_hist(:, t_idx) = v_sat;

    trq_hist(t_idx) = trqCalc(i);

    L1_hist(:, t_idx) = nnOpt1.Lambda;
    th1_hist(:, t_idx) = nn_V_norm_cal(nn1.V, nnOpt1);
    L2_hist(:, t_idx) = nnOpt2.Lambda;
    th2_hist(:, t_idx) = nn_V_norm_cal(nn2.V, nnOpt2);

    if mod(t_idx, rpt_dt/dt) == 0
        % fprintf('t_idx: %d\n', t_idx)
        fprintf("Time : %.5f\n", t(t_idx));
    end 
    
end

%%
fig_height = 210; fig_width = 450;

figure(1); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];
plot(t, wr_hist, 'b', 'LineWidth', 2, 'LineStyle', '-'); hold on;
plot(t, w_hist, 'r', 'LineWidth', 2, 'LineStyle','-.'); hold on;
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Speed [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
maxVal = max(wr_hist); minVal = min(wr_hist); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);

figure(2); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(2,1);

nexttile;
plot(t, ir_hist(1,:), 'b', 'LineWidth', 2, 'LineStyle', '-'); hold on;
plot(t, i_hist(1,:), 'r', 'LineWidth', 2, 'LineStyle', '-.'); hold on;
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('id [A]', 'FontSize', 12, 'Interpreter', 'latex');
maxVal = max(ir_hist(1,:)); minVal = min(ir_hist(1,:)); 
% maxVal = nnOpt1.cstr.u_ball; minVal = -nnOpt1.cstr.u_ball;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);

nexttile
plot(t, ir_hist(2,:), 'b', 'LineWidth', 2, 'LineStyle', '-'); hold on;
plot(t, i_hist(2,:), 'r', 'LineWidth', 2, 'LineStyle', '-.'); hold on;
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('iq [A]', 'FontSize', 12, 'Interpreter', 'latex');
maxVal = max(ir_hist(2,:)); minVal = min(ir_hist(2,:)); 
% maxVal = nnOpt1.cstr.u_ball; minVal = -nnOpt1.cstr.u_ball;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);

figure(3); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];
tiledlayout(1,2);

nexttile
norm_i = sqrt(sum(i_hist.^2)); norm_id = sqrt(sum(ir_hist.^2));
plot(t, ones(size(t))*nnOpt1.cstr.u_ball, 'k', 'LineWidth',2); hold on
plot(t, norm_id, 'b', 'LineWidth', 2, 'LineStyle', '-'); hold on
plot(t, norm_i, 'r', 'LineWidth', 2, 'LineStyle', '-.'); hold on
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\Vert i\Vert$ [A]', 'FontSize', 12, 'Interpreter', 'latex');

nexttile
ang = 0:0.01:2*pi;
plot(nnOpt1.cstr.u_ball*cos(ang), nnOpt1.cstr.u_ball*sin(ang), 'k', 'LineWidth', 2); hold on;
plot(i_hist(1,:), i_hist(2,:), 'r', 'LineWidth', 2); hold on;
plot(ir_hist(1,:), ir_hist(2,:), 'b', 'LineWidth', 2); hold on;
grid on;
xlabel('id [A]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('iq [A]', 'FontSize', 12, 'Interpreter', 'latex');

figure(4); clf;
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(2,1);

nexttile;
plot(t, v_hist(1,:), 'r', 'LineWidth', 2); hold on;
plot(t, vSat_hist(1,:), 'b', 'LineWidth', 2); hold on;
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('vd [V]', 'FontSize', 12, 'Interpreter', 'latex');
maxVal = max(vSat_hist(1,:)); minVal = min(vSat_hist(1,:)); 
% maxVal = nnOpt2.cstr.u_ball; minVal = -nnOpt2.cstr.u_ball;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);

nexttile
plot(t, v_hist(2,:), 'r', 'LineWidth', 2); hold on;
plot(t, vSat_hist(2,:), 'b', 'LineWidth', 2); hold on;
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('vq [V]', 'FontSize', 12, 'Interpreter', 'latex');
maxVal = max(vSat_hist(2,:)); minVal = min(vSat_hist(2,:)); 
% maxVal = nnOpt2.cstr.u_ball; minVal = -nnOpt2.cstr.u_ball;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);

figure(5);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(1,2);

nexttile
norm_v = sqrt(sum(v_hist.^2));
plot(t, ones(size(t))*nnOpt2.cstr.u_ball, 'k', 'LineWidth',2); hold on
plot(t, norm_v, 'r', 'LineWidth', 2, 'LineStyle', '-'); hold on
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\Vert v\Vert$ [V]', 'FontSize', 12, 'Interpreter', 'latex');

nexttile
ang = 0:0.01:2*pi;
plot(nnOpt2.cstr.u_ball*cos(ang), nnOpt2.cstr.u_ball*sin(ang), 'k', 'LineWidth', 2); hold on;
plot(v_hist(1,:), v_hist(2,:), 'r', 'LineWidth', 2); hold on;
grid on
xlabel('vd [V]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('vq [V]', 'FontSize', 12, 'Interpreter', 'latex');

figure(6);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(1,2);

nexttile
plot(t, th1_hist(1,:), 'r', 'LineWidth', 2); hold on;
plot(t, th1_hist(2,:), 'b', 'LineWidth', 2); hold on;
grid on;

nexttile
semilogy(t, L1_hist(1,:), 'r', 'LineWidth', 2); hold on;
semilogy(t, L1_hist(2,:), 'g', 'LineWidth', 2); hold on;
semilogy(t, L1_hist(3,:), 'b', 'LineWidth', 2); hold on;
semilogy(t, L1_hist(3,:), 'LineWidth', 2); hold on;
grid on;

figure(7);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];
tl = tiledlayout(1,2);

nexttile
plot(t, th2_hist(1,:), 'r', 'LineWidth', 2); hold on;
plot(t, th2_hist(2,:), 'b', 'LineWidth', 2); hold on;
grid on;

nexttile
semilogy(t, L2_hist(1,:), 'r', 'LineWidth', 2); hold on;
semilogy(t, L2_hist(2,:), 'g', 'LineWidth', 2); hold on;
semilogy(t, L2_hist(3,:), 'b', 'LineWidth', 2); hold on;
semilogy(t, L2_hist(4,:), 'LineWidth', 2); hold on;
semilogy(t, L2_hist(5,:), 'LineWidth', 2); hold on;
grid on;

figure(8);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];
plot(t, -ones(size(t))*nnOpt2.cstr.tau_ball, 'k', 'LineWidth',2); hold on
plot(t, trq_hist, 'r', 'LineWidth', 2);
grid on;
xlabel('Time [s]', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('$\tau$ [Nm]', 'FontSize', 12, 'Interpreter', 'latex');

%% LOCAL FUNCTION
function grad = sysGrad_1st(w, i)
    J = [0 -1;1 0];

    L = [
        0.45     0;
        0     0.66
    ] * 1e-3;
    R = [
        25     0;
        0     25;
    ] * 1e-3;

    np = 8;             
    kappa = 2/3;          
    Theta = 0.005;   

    bias_psi = [0.0563; 0];
    % 
    % k1 = 1.5*np*bias_psi(1);
    % k2 = 1.5*np*(L(1,1) - L(2,2));
    % k3 = k1/(2*k2);

    psi = L*i + bias_psi;

    trq = (2*np)/(3*kappa^2) * i'*J*psi;
    % trq = (2*np)/(3*kappa^2) * i'*J*psi;
    % trq = (k1+k2*i(1)) * i(2);
    trq_l = 0;
    grad = 1/Theta * (trq - trq_l);   


end


function grad = sysGrad_linear(w, i, u)
    J = [0 -1;1 0];

    L = [
        0.45     0;
        0     0.66
    ] * 1e-3;
    R = [
        25     0;
        0     25;
    ] * 1e-3;

    np = 8;               
    kappa = 2/3;          
    Theta = 0.005;   

    bias_psi = [0.0563; 0];

    k1 = 1.5*np*bias_psi(1);
    k2 = 1.5*np*(L(1,1) - L(2,2));
    k3 = k1/(2*k2);

    psi = L*i + bias_psi;

    trq = (2*np)/(3*kappa^2) * i'*J*psi;
    % trq = (k1+k2*i(1)) * i(2);
    % trq_l = 1*tanh(w);
    trq_l = 0;
    grad1 = 1/Theta * (trq - trq_l);   

    inv_L = matInv22(L);
    grad2 = inv_L * (-R*i - w/np*J*psi + u);

    grad = [grad1; grad2];

end

function grad = sysGrad_lookup(w, i, u, param)
    J = [0 -1; 1 0];

    % current vectors / grids
    isd = param.machine.psi.s.arg.isd;
    isq = param.machine.psi.s.arg.isq;

    % flux linkages
    PSISD = param.machine.psi.s.d;
    PSISQ = param.machine.psi.s.q;

    % differential inductances
    LSDD = param.machine.L.s.dd;
    LSDQ = param.machine.L.s.dq;
    LSQD = param.machine.L.s.qd;
    LSQQ = param.machine.L.s.qq;

    % dummy data
    R = eye(2) * param.machine.Rs;
    
    np = param.machine.nP;               
    kappa = param.machine.kappa;          
    Theta = param.machine.ThetaM;

    L_dd = interpolate_from_table(i(1), i(2), isd, isq, LSDD);
    L_dq = interpolate_from_table(i(1), i(2), isd, isq, LSDQ);
    L_qd = interpolate_from_table(i(1), i(2), isd, isq, LSQD);
    L_qq = interpolate_from_table(i(1), i(2), isd, isq, LSQQ);
    L = [L_dd, L_dq; L_qd, L_qq];
    L = L*1;

    psi_d = interpolate_from_table(i(1), i(2), isd, isq, PSISD);
    psi_q = interpolate_from_table(i(1), i(2), isd, isq, PSISQ);
    psi = [psi_d; psi_q];
    psi = psi*1;


    trq = (2*np)/(3*kappa^2) * i'*J*psi;
    % trq = (k1+k2*i(1)) * i(2);
    % trq_l = 1*sign(w);
    trq_l = 0;
    grad1 = 1/Theta * (trq - trq_l);   

    inv_L = matInv22(L);
    grad2 = inv_L * (-R*i - w/np*J*psi + u);

    grad = [grad1; grad2];
end

function trq = trqCalc(i)

    L = [
        0.45     0;
        0     0.66
    ] * 1e-3;
    np = 8;               
    bias_psi = [0.0563; 0];

    k1 = 1.5*np*bias_psi(1);
    k2 = 1.5*np*(L(1,1) - L(2,2));
    k3 = k1/(2*k2);

    trq = (k1+k2*i(1)) * i(2);
end

function inv_M = matInv22(M)
    det = M(1,1)*M(2,2) - M(1,2)*M(2,1);
    assert(det ~= 0, 'Matrix is singular and cannot be inverted')

    inv_M = (1/det) * [
        +M(2,2), -M(1,2); 
        -M(2,1), +M(1,1)
    ];
end

function y = linear_interpolate(x, x1, x2, y1, y2)
    % y = (y2-y1)/(x2-x1) * (x - x1) + y1     
    y = (y2-y1)/(x2-x1) * (x-x1) + y1;
end

function z = bilinear_interpolate(x,y, x1,x2,y1,y2, z11,z12,z21,z22)
    z1 = linear_interpolate(x, x1, x2, z12, z22);
    z2 = linear_interpolate(x, x1, x2, z11, z21);
    z = linear_interpolate(y, y1, y2, z2, z1);
end

function y = interpolate_from_table(x,y, X,Y, Z)
    x_idx = find(X >= x); x_idx = x_idx(1);
    y_idx = find(Y >= y); y_idx = y_idx(1);

    x_upper = X(x_idx); x_lower = X(x_idx-1);
    y_upper = Y(y_idx); y_lower = Y(y_idx-1);

    z_ul = Z(x_idx-1, y_idx-1);
    z_uu = Z(x_idx, y_idx-1);
    z_ll = Z(x_idx-1, y_idx);
    z_lu = Z(x_idx, y_idx);

    y = bilinear_interpolate(x, y, x_lower, x_upper, y_lower, y_upper, z_ul, z_uu, z_ll, z_lu);
end