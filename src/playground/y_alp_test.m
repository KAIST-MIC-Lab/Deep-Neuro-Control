figure(1); clf;

B = eye(2);
inv_R = eye(2);
SDC = [
    -10.0000   1.0000
    1   -10.0000
];
e = [3;-5]*1e0;
x_num = 2;
dt = 1e-2;
M_pre = eye(2);

alpha = 1e1;

% contraction constraint
H = -2*B*inv_R*B';
f = transpose(e'*(1/dt*eye(x_num)+2*SDC'+2*alpha*eye(x_num)));
g = -e'*M_pre*e/dt;

J_max = 1/4*f'*inv(H)*f + g;
fprintf("Maximum of the contraction constraint: %.2f\n", J_max);

% H = -H; f = -f; g = -g;
% fp = fimplicit(@(x,y) H(1,1)*x.^2 + (H(1,2)+H(2,1))*x.*y + H(2,2)*y.^2 + f(1)*x + f(2)*y + g, 'Color', 'black', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'blue')

func = @(x,y) H(1,1)*x.^2 + (H(1,2)+H(2,1))*x.*y + H(2,2)*y.^2 + f(1)*x + f(2)*y + g;
% func = @(x,y) [x;y]' * H * [x;y] + f'*[x;y] + g;

max_X = 1000; dx = 0.5;
[X,Y] = meshgrid(-max_X:dx:max_X, -max_X:dx:max_X);
Z = func(X,Y);
% plot3(x,y, func(x',y), 'k'); hold on
contour3(X,Y, Z, 'LineColor', 'black', 'LineWidth', 1);

grid on

xlim([-max_X max_X])
ylim([-max_X max_X])
zlim([-500 500])

