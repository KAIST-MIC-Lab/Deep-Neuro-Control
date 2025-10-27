clear
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));
addpath(genpath([pwd filesep 'mosek/11.0/tools/platform/osxaarch64/bin']));
addpath(genpath([pwd filesep 'mosek/11.0/toolbox/r2022bom']));
% C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin

%%
SEED = 1;
rng(SEED);

%% 
x_num = 2;
dt = 1e-3;

%% GIVEN VALUES
% SDC = diag([-10,-1]);
B = diag([1,3])*1e0;
inv_R = eye(2);

SDC = [
    -10.0000   1.0000
    5   -3.0000
];

% W_pre = eye(x_num)*1e2;
% W_pre = (rand(x_num)+diag([2;3])) * 1e0; W_pre = W_pre'+W_pre;
M_pre = [
    15  6
    6 10
]*1e0;
W_pre = inv(M_pre);
mu_pre = 1;
% M_pre = inv(W_pre)*mu_pre;

e = [3;-5]*1e0;
alpha = 1e0;

%% 
ud = [1;-.5];

m_max = 1e1;
m_min = 1e-3;
u_max = 3.5;

if ud'*ud - u_max^2 > 0
    error("desired input exceeds the limit")
end

%% VISUALIZATION (y-space)
figure(1);clf; 

% contraction constraint
H = -2*B*inv_R*B';
f = transpose(e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num)));
g = -e'*M_pre*e/dt;
% H = -H; f = -f; g = -g;
fp = fimplicit(@(x,y) H(1,1)*x.^2 + (H(1,2)+H(2,1))*x.*y + H(2,2)*y.^2 + f(1)*x + f(2)*y + g, 'Color', 'black', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'blue')

% input saturation constraint
H = B*inv_R*inv_R*B';
f = transpose(-2*ud'*inv_R*B');
g = ud'*ud - u_max^2;
fp = fimplicit(@(x,y) H(1,1)*x.^2 + (H(1,2)+H(2,1))*x.*y + H(2,2)*y.^2 + f(1)*x + f(2)*y + g, 'Color', 'magenta', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'red')

% I points
% mu_list = m_min:1e-1:m_max;
% for mu_idx = 1:1:length(mu_list)
%     mu = mu_list(mu_idx);
%     y_I = eye(x_num)*mu*e;
% 
%     u_I = ud - inv_R*B'*y_I;
% 
%     plot([0, y_I(1)], [0, y_I(2)], 'Color', [0.5 0.5 0.5], 'LineStyle','--', 'LineWidth', 1, 'Marker','o'); hold on;
%     % text(y_I(1), y_I(2), "I: u* = ["+num2str(u_I(1), '%.2f')+", "+num2str(u_I(2), '%.2f')+"]", 'Color', [0.5 0.5 0.5]);
% end

daspect([1 1 1])
grid on

%% EXISTING PROBLEM
W_bar = sdpvar(x_num ,x_num); assign(W_bar, W_pre)
mu = sdpvar(1,1); assign(mu, mu_pre);
% X = sdpvar(1,1); assign(X, mu_pre/m_min);
X = sdpvar(1,1); assign(X, 100);

lbd = 1e-50;

% mu = 1;
con = [ 
    -(W_bar-W_pre)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
    eye(x_num) <= W_bar;
    W_bar <= X * eye(x_num);
    0 <= mu
    % X == mu/m_min;
    % m_min <= mu;
    % 1e-99 <= mu;
    1 <= X;
    % X == mu/1e-99;
];

obj = X + lbd * mu;

ops = sdpsettings('solver','sedumi','verbose',0);
sol = optimize(con, obj, ops);

if sol.problem ~= 0
    disp('Something went wrong!');
    error("Optimization failed: "+yalmiperror(sol.problem));
    return;
end

W_bar = value(W_bar);
mu = value(mu);
X = value(X);
M = mu * inv(W_bar);

y_existing = M*e;
u_existing = ud - inv_R*B'*y_existing;

plot([0, y_existing(1)], [0, y_existing(2)], 'Color', 'red', 'LineStyle','-', 'LineWidth', 2, 'Marker','o');
text(y_existing(1), y_existing(2), "Existing "+newline+"u* = ["+num2str(u_existing(1), '%.2f')+", "+num2str(u_existing(2), '%.2f')+"]", 'Color', 'red');

fprintf("EXISTING PROBLEM SOLUTION: u* = ["+num2str(u_existing(1), '%.2f')+", "+num2str(u_existing(2), '%.2f')+"\n");
fprintf("W_bar = \n"); disp(W_bar);
fprintf("M = \n"); disp(M);
fprintf("mu = %.4f, X = %.4f\n", mu, X);

eig(-(W_bar-W_pre)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' + 2*alpha*W_bar)
eig((M-M_pre)/dt + (SDC'*M+M*SDC) - 2*M*B*inv_R*B'*M + 2*alpha*M)

%% PROPOSED
M = sdpvar(x_num ,x_num); assign(M, M_pre);
y = sdpvar(2,1);

con = [
    -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(2)+2*SDC'+2*alpha*eye(2))*y - (e'*M_pre*e)/dt <= 0;
    % y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
    M*e == y;
];

obj = norm(M-eye(2));

ops = sdpsettings('solver', 'fmincon', 'verbose', 0, 'fmincon.Algorithm', 'sqp');
sol = optimize(con, obj, ops);
if sol.problem ~= 0
    disp('Something went wrong!');
    error("Optimization failed: "+yalmiperror(sol.problem));
    return;
end

M = value(M);
y_proposed = M*e;
u_proposed = ud - inv_R*B'*y_proposed;

fprintf("PROPOSED PROBLEM SOLUTION (without constraint):\n");
fprintf("M = \n"); disp(M);
% fprintf("y_proposed = \n"); disp(y_proposed);
plot([0, y_proposed(1)], [0, y_proposed(2)], 'Color', 'cyan', 'LineStyle','-', 'LineWidth', 2, 'Marker','o');

text(y_proposed(1), y_proposed(2), "Proposed (without cstr)"+newline+"u* = ["+num2str(u_proposed(1), '%.2f')+", "+num2str(u_proposed(2), '%.2f')+"]", 'Color', 'cyan');

%% PROPOSED (with constraint)
M = sdpvar(x_num ,x_num); assign(M, M_pre);
y = sdpvar(2,1);

con = [
    -y'*(2*B*inv_R*B')*y + e'*(1/dt*eye(2)+2*SDC'+2*alpha*eye(2))*y - (e'*M_pre*e)/dt <= 0;
    y'*(B*inv_R*inv_R*B')*y + (-2*ud'*inv_R*B')*y + (ud'*ud - u_max^2) <= 0;
    M*e == y;
];

obj = norm(M-eye(2));

ops = sdpsettings('solver', 'fmincon', 'verbose', 0, 'fmincon.Algorithm', 'sqp');
sol = optimize(con, obj, ops);
if sol.problem ~= 0
    disp('Something went wrong!');
    error("Optimization failed: "+yalmiperror(sol.problem));
    return;
end

M = value(M);
y_proposed_cstr = M*e;
u_proposed_cstr = ud - inv_R*B'*y_proposed;

fprintf("PROPOSED PROBLEM SOLUTION (with constraint):\n");
fprintf("M = \n"); disp(M);
fprintf("")
% fprintf("y_proposed_cstr = \n"); disp(y_proposed_cstr);
plot([0, y_proposed_cstr(1)], [0, y_proposed_cstr(2)], 'Color', 'blue', 'LineStyle','-', 'LineWidth', 2, 'Marker','o');

text(y_proposed_cstr(1), y_proposed_cstr(2), "Proposed (with cstr) "+newline+"u* = ["+num2str(u_proposed_cstr(1), '%.2f')+", "+num2str(u_proposed_cstr(2), '%.2f')+"]", 'Color', 'blue');
