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
B = diag([1,2])*1e0;
inv_R = eye(2);

SDC = [
    -10.0000   7.0000
    5   -1.0000
];

% W_pre = eye(x_num)*1e2;
W_pre = (rand(x_num)+diag([2;3])) * 1e0; W_pre = W_pre'+W_pre;
mu_pre = 1;
M_pre = inv(W_pre)*mu_pre;

e = [3;5]*1e1;
alpha = 5e1;

%% 
ud = [1;0];

m_min = 1e-3;
u_max = 35;

if ud'*ud - u_max^2 > 0
    error("desired input exceeds the limit")
end

%% VISUALIZATION (y-space)
figure(1);clf; 

% contraction constraint
H = -2*B*inv_R*B';
f = transpose(e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num)));
g = -e'*M_pre*e/dt;
fp = fimplicit(@(x,y) [x;y]'*H*[x;y]+f'*[x;y]+g, 'Color', 'black', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'blue')

% input saturation constraint
H = B*inv_R*inv_R*B';
f = transpose(-2*ud'*inv_R*B');
g = ud'*ud - u_max^2;
fp = fimplicit(@(x,y) [x;y]'*H*[x;y]+f'*[x;y]+g, 'Color', 'magenta', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'red')

daspect([1 1 1])
grid on

%% EXISTING PROBLEM
W_bar = sdpvar(x_num ,x_num); assign(W_bar, W_pre)
mu = sdpvar(1,1); assign(mu, mu_pre);
X = sdpvar(1,1); assign(X, mu_pre/m_min);

lbd = 5e-2;

con = [ 
    -(W_bar-W_pre)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
    eye(x_num) <= W_bar;
    W_bar <= X * eye(x_num);
    0 <= mu
    % X == mu/m_min;
    % m_min <= mu;
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

y_existing = mu*inv(W_bar)*e;
u_existing = ud - inv_R*B'*y_existing;

plot([0, y_existing(1)], [0, y_existing(2)], 'Color', 'red', 'LineStyle','-', 'LineWidth', 2, 'Marker','o');
text(y_existing(1), y_existing(2), "Existing "+newline+"u* = ["+num2str(u_existing(1), '%.2f')+", "+num2str(u_existing(2), '%.2f')+"]", 'Color', 'red');

fprintf("EXISTING PROBLEM SOLUTION: u* = ["+num2str(u_existing(1), '%.2f')+", "+num2str(u_existing(2), '%.2f')+"\n");
fprintf("W_bar = \n"); disp(W_bar);
fprintf("M = \n"); disp(mu*inv(W_bar));
fprintf("mu = %.4f, X = %.4f\n", mu, X);
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
