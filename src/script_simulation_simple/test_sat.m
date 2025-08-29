clear

A = diag([1e-1 1e-1]);
B = [5 1; 2 3];

e = [3;-4];

ud = [1.5;0.5];

%%

u_test = ud - A*B'*e;
fprintf("test u is %f\n", norm(u_test))

max_u = 1e-2;

%%
% u = A*B'* M *e;
AUX = sdpvar(2, 2);
W_bar = sdpvar(2, 2);

con = [
    [max_u^2                 (W_bar*inv(A*B')*ud - e)'
    W_bar*inv(A*B')*ud - e      AUX] >= 0
];

% con = [con,
%     [AUX (W_bar*inv(A*B'))'
%     (W_bar*inv(A*B')) eye(2)] >= 0,
% ];

con = [con,
    [AUX (W_bar)'
    (W_bar) inv(A*B')*inv(A*B')'] >= 0,
];

%%
ops = sdpsettings('verbose', 0);
ops = sdpsettings(ops, 'debug',0);


ops = sdpsettings(ops,'solver', 'sedumi');

sol = optimize(con, 1, ops);

%%
W_bar = value(W_bar);
AUX = value(AUX);
fprintf("LMI result: %d\n", sol.problem)

u_test = ud - A*B'*inv(W_bar)*e;
fprintf("test u is %f\n", norm(u_test))
