clear

%% 
x_num = 2;
dt = 1e-3;

%% GIVEN VALUES
% SDC = -diag([1,2]);
B = diag([1,1])*1e0;
inv_R = eye(2);

SDC = [
  -10.0000   7.0000
    5   -1.0000
   ];

% W_pre = eye(x_num)*1e-1;
W_pre = (rand(x_num)+diag([2;3])) * 1e0; W_pre = W_pre'+W_pre;
mu_pre = 10;
M_pre = inv(W_pre)*mu_pre;

e = rand(2,1)*-2;
% e = [1;1];
alpha = 1e0;

%% 
ud = [1;0.5];

m_min = 1e-5;
u_max = 3.5;

%% VISUALIZATION
% u-space; u = ud - inv_R*transpose(B)*M*e; M=mu*inv(W_bar);

angles = 0:0.01:pi/2;

figure(1);clf; axis equal;

for ang = angles
    r = sdpvar(1,1); assign(r, 0.001);
    P = sdpvar(2,1);
    M = sdpvar(2,2);

    con =  [
        r >= 0;
        P(1) == r*cos(ang);
        P(2) == r*sin(ang);
        P == ud - inv_R*B'*M*e;
        (M-M_pre)/dt + (M*SDC'+SDC*M) - 2*M*B*inv_R*B'*M <= -2*alpha*M;
        m_min*eye(2) <= M;
        ];

    obj = -r;

    % ops = sdpsettings('solver','sedumi','verbose',0);
    ops = sdpsettings('solver', 'fmincon', 'verbose', 0, 'fmincon.Algorithm', 'sqp');

    sol = optimize(con, obj, ops);

    if sol.problem ~= 0
        % disp('Something went wrong!');
        % fprintf("Optimization failed: "+yalmiperror(sol.problem));

    else
        r_opt = value(r);
        P_opt = value(P);
        M_opt = value(M);
    
        plot(P_opt(1), P_opt(2), 'Color', 'blue', 'LineStyle', 'none', 'Marker','o');hold on;
    end

    plot(u_max*cos(ang), u_max*sin(ang), 'Color', 'magenta', 'LineStyle', 'none', 'Marker','o');
end

grid on;

% ud
plot([0,ud(1)],[0,ud(2)], 'Color','red');
text(ud(1), ud(2), "ud = ["+num2str(ud(1), '%.2f')+", "+num2str(ud(2), '%.2f')+"]", 'Color','red');
%

% figure(1); clf;
% theta = 0:0.01:2*pi;
% x = cos(theta);

%%

%% EXISTING PROBLEM
W_bar = sdpvar(x_num ,x_num); assign(W_bar, W_pre)
mu = sdpvar(1,1); assign(mu, mu_pre);
X = sdpvar(1,1); assign(X, mu_pre/m_min);

lbd = 1e-3;

con = [ 
    -(W_bar-W_pre)/dt + (W_bar*SDC' + SDC*W_bar) - 2*mu*B*inv_R*B' <= -2*alpha*W_bar;
    eye(x_num) <= W_bar;
    W_bar <= X * eye(x_num);
    % 0 <= mu
    X == mu/m_min;
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

fprintf("EXISTING PROBLEM SOLUTION:\n");
fprintf("W_bar = \n"); disp(W_bar);
fprintf("mu = %.4f\n", mu);
fprintf("X = %.4f\n", X);

ud_existing = ud - inv_R*B'*mu*inv(W_bar)*e;
    
plot([ud(1), ud_existing(1)], [ud(2), ud_existing(2)], 'Color', 'green', 'LineStyle','--');
text(ud_existing(1), ud_existing(2), "Existing "+newline+"u* = ["+num2str(ud_existing(1), '%.2f')+", "+num2str(ud_existing(2), '%.2f')+"]");

%% PROPOSED
M = sdpvar(x_num ,x_num); assign(M, M_pre);
y = sdpvar(2,1);

con = [
    -y'*(B*inv_R*B')*y + e'*(1/dt*eye(2)-2*SDC'+2*alpha*eye(2))*y - (e'*M_pre*e)/dt <= 0;
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
ud_proposed = ud - inv_R*B'*M*e;

fprintf("PROPOSED PROBLEM SOLUTION:\n");
fprintf("M = \n"); disp(M);
fprintf("ud_proposed = \n"); disp(ud_proposed);
plot([ud(1), ud_proposed(1)], [ud(2), ud_proposed(2)], 'Color', 'cyan', 'LineStyle','--');

text(ud_proposed(1), ud_proposed(2), "Proposed "+newline+"u* = ["+num2str(ud_proposed(1), '%.2f')+", "+num2str(ud_proposed(2), '%.2f')+"]");

%% PROPOSED (with constraint)
M = sdpvar(x_num ,x_num); assign(M, M_pre);
y = sdpvar(2,1);

con = [
    -y'*(B*inv_R*B')*y + e'*(1/dt*eye(2)-2*SDC'+2*alpha*eye(2))*y - (e'*M_pre*e)/dt <= 0;
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
ud_proposed = ud - inv_R*B'*M*e;

fprintf("PROPOSED PROBLEM SOLUTION:\n");
fprintf("M = \n"); disp(M);
fprintf("ud_proposed = \n"); disp(ud_proposed);
plot([ud(1), ud_proposed(1)], [ud(2), ud_proposed(2)], 'Color', 'cyan', 'LineStyle','--');

text(ud_proposed(1), ud_proposed(2), "Proposed (with cstr) "+newline+"u* = ["+num2str(ud_proposed(1), '%.2f')+", "+num2str(ud_proposed(2), '%.2f')+"]");
