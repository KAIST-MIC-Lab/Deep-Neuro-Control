clear

SAVE_FLAG = 1;

%% 
MAT_NAME = "1-09-2025_17-01-19";
load("results/"+MAT_NAME+".mat");

%%
param.max_u = 1.5e2;
num_x = 3; num_u = 3;
num_sample = length(sys);
T = t(end);

%% FIGURE SETTING
POSITION_FLAG = 1; % it will plot fiugures in the same position

font_size = 16;
line_width = 2;
lgd_size = 16;
fig_height = 200; 
fig_width = 450;

% For papers
% font_size = 32;
% line_width = 2;
% lgd_size = 28;
% fig_height = 300; 
% fig_width = 800;

%% COLORS
c_list = [
    0 0 1       % BLUE
    0 1 1       % CYAN
    .5 .5 .5    % GRAY
];

%% MAIN PLOT
% ============================================
%     figure 1: STATE VARIABLE
% ============================================
figure(1); clf; 
hF = gcf; 
hF.Position(3:4) = [1600, 200];
tl = tiledlayout(1, 7, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     STATE VARIABLE
% ============================================
for x_idx = 1:1:num_x
    nexttile([1,2]);

    maxVal = 0; minVal = 0;
    for idx = 1:1:num_sample
        c = c_list(idx, :);

        plot(t, sys{idx}.x_hist(x_idx,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
        plot(t, sys{idx}.x_non_hist(x_idx,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on

        maxVal = max(maxVal, max([sys{idx}.x_hist(x_idx,:), sys{idx}.x_non_hist(x_idx,:), xd_hist(x_idx,:)]));
        minVal = min(minVal, min([sys{idx}.x_hist(x_idx,:), sys{idx}.x_non_hist(x_idx,:), xd_hist(x_idx,:)]));
    end
    plot(t, xd_hist(x_idx,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

    grid on; grid minor;
    xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
    ylabel("$x_"+x_idx+"$", 'FontSize', font_size, 'Interpreter', 'latex');
    len = maxVal-minVal; ratio = .1;
    if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
    xlim([0 T])

    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';
end


% ============================================
%     STATE VARIABLE (Bird-eye view)
% ============================================
nexttile;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    plot3(sys{idx}.x_hist(1,:), sys{idx}.x_hist(2,:), sys{idx}.x_hist(3,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    plot3(sys{idx}.x_non_hist(1,:), sys{idx}.x_non_hist(2,:), sys{idx}.x_non_hist(3,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on
end
plot3(xd_hist(1,:), xd_hist(2,:), xd_hist(3,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

grid on; grid minor;
xlabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
zlabel('$x_3$', 'FontSize', font_size, 'Interpreter', 'latex');

max_x = 150;
xlim([-max_x max_x]);
ylim([-max_x max_x]);
zlim([-max_x max_x]);

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';
% pbaspect([1 1 1])
% axis equal
% axis square

% ============================================
%     figure 2: CONTROL INPUT
% ============================================
figure(2); clf; 
hF = gcf; 
hF.Position(3:4) = [1600, 200];
tl = tiledlayout(1, 7, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     CONTROL INPUT
% ============================================
for u_idx = 1:1:num_u
    nexttile([1,2]);

    maxVal = 0; minVal = 0;
    for idx = 1:1:num_sample
        c = c_list(idx, :);

        plot(t, sys{idx}.u_hist(u_idx,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
        plot(t, sys{idx}.uSat_hist(u_idx,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on

        maxVal = max(maxVal, max([sys{idx}.u_hist(u_idx,:), sys{idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
        minVal = min(minVal, min([sys{idx}.u_hist(u_idx,:), sys{idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
    end
    plot(t, ud_hist(u_idx,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

    grid on; grid minor;
    xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
    ylabel("$u_"+u_idx+"$", 'FontSize', font_size, 'Interpreter', 'latex');
    len = maxVal-minVal; ratio = .1;
    if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
    xlim([0 T])

    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';
end

% ============================================
%     CONTROL INPUT LOCI
% ============================================
nexttile;
max_u = param.max_u;

for idx = 1:1:num_sample
    c = c_list(idx, :);

    plot3(sys{idx}.u_hist(1,:), sys{idx}.u_hist(2,:), sys{idx}.u_hist(3,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    plot3(sys{idx}.uSat_hist(1,:), sys{idx}.uSat_hist(2,:), sys{idx}.uSat_hist(3,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on
end
[x y z] = sphere(128);
h = surfl(max_u*x, max_u*y, max_u*z); 
set(h, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% shading interp
% plot3(ud_hist(1,:), ud_hist(2,:), ud_hist(3,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

% rad = -pi:0.1:pi;
% plot3(max_u*sin(rad), max_u*cos(rad), zeros(size(rad)), "Color", "black")

grid on; grid minor;
xlabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
zlabel('$u_3$', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0
%     ylim([minVal-len*ratio maxVal+len*ratio]);
%     xlim([minVal-len*ratio maxVal+len*ratio]);
% end

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';
% pbaspect([1 1 1])
axis equal
% axis square

%% 
% ============================================
%     figure 3: CONTROL INPUT NORM
% ============================================
figure(3); clf; 
hF = gcf; 
hF.Position(3:4) = [400, 200];
% tl = tiledlayout(1, 4, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     Control norm
% ============================================
% nexttile;

plot([-10, T+10], param.max_u * [1 1], "Color", "black", "LineWidth", line_width, "LineStyle", "-"); hold on
maxVal = 0; minVal = 0;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    u_norm = time_norm(sys{idx}.u_hist);
    plot(t, u_norm, "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(u_norm));
    minVal = min(minVal, min(u_norm));
end

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel("$\Vert u\Vert$", 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     figure 4: TRACKING ERROR
% ============================================
figure(4); clf; 
hF = gcf; 
hF.Position(3:4) = [400, 200];
tl = tiledlayout(1, 4, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     Tracking error 
% ============================================
maxVal = 0; minVal = 0;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    err = time_norm(sys{idx}.x_hist-xd_hist);
    plot(t, err, "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(err));
    minVal = min(minVal, min(err));
end

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel("$e$", 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';


% % ============================================
% %     X
% % ============================================
% nexttile;

% maxVal = 0; minVal = 0;
% for idx = 1:1:num_sample
%     c = c_list(idx, :);

%     plot(t, sys{idx}.M_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    
%     maxVal = max(maxVal, max(sys{idx}.M_hist(1,:)));
%     minVal = min(minVal, min(sys{idx}.M_hist(1,:)));
% end

% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$X$', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
% xlim([0 T])

% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';

% % ============================================
% %     optDone
% % ============================================
% nexttile;

% maxVal = 0; minVal = 0;
% for idx = 1:1:num_sample
%     c = c_list(idx, :);

%     plot(t, sys{idx}.optDone_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

%     maxVal = max(maxVal, max(sys{idx}.optDone_hist(1,:)));
%     minVal = min(minVal, min(sys{idx}.optDone_hist(1,:)));
% end

% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('Done', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
% xlim([0 T])

% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';

% % ============================================
% %     mu
% % ============================================
% nexttile;

% maxVal = 0; minVal = 0;
% for idx = 1:1:num_sample
%     c = c_list(idx, :);

%     plot(t, sys{idx}.mu_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

%     maxVal = max(maxVal, max(sys{idx}.mu_hist(1,:)));
%     minVal = min(minVal, min(sys{idx}.mu_hist(1,:)));
% end

% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('mu', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
% xlim([0 T])

% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';



%% SAVE FIGURES
if SAVE_FLAG
    [~,~] = mkdir("figures/compare");

    for idx = 1:1:10

        f_name = "figures/compare/Fig" + string(idx);

        saveas(figure(idx), f_name + ".png")

        figure(idx);
        % set(gcf, 'Position', [0, 0, fig_width, fig_height]); % [left, bottom, width, height] 
        exportgraphics(gcf, f_name+'.eps', 'ContentType', 'vector')
        % exportgraphics(figure(idx), f_name+'.eps',"Padding","figure")
    
        % matlab2tikz(char(f_name+".tex"))
    end
end

%% NUMERICAL ANALYSIS
% ctrl_dt = 1/250;
% sim_dt = ctrl_dt / 1000;

% e1 = [
%     transpose(data1.q1.Data-data1.r1.Data);
%     transpose(data1.q2.Data-data1.r2.Data)
% ];
% e1 = e1(:, (data1.r1.Time >= idleTime1 & data1.r1.Time <= idleTime1+T));
% e2 = [
%     transpose(data2.q1.Data-data2.r1.Data);
%     transpose(data2.q2.Data-data2.r2.Data)
% ];
% e2 = e2(:, (data2.r1.Time >= idleTime1 & data2.r1.Time <= idleTime1+T));
% e3 = [
%     transpose(data3.q1.Data-data3.r1.Data);
%     transpose(data3.q2.Data-data3.r2.Data)
% ];
% e3 = e3(:, (data3.r1.Time >= idleTime1 & data3.r1.Time <= idleTime1+T));
% e4 = [
%     transpose(data4.q1.Data-data4.r1.Data);
%     transpose(data4.q2.Data-data4.r2.Data)
% ];
% e4 = e4(:, (data4.r1.Time >= idleTime1 & data4.r1.Time <= idleTime1+T));

% ep_idx = floor(size(e1,2)/2);

% e11_ep1 = e1(1,1:ep_idx);
% e12_ep1 = e1(2,1:ep_idx);
% e11_ep2 = e1(1,ep_idx+1:end);
% e12_ep2 = e1(2,ep_idx+1:end);
% e21_ep1 = e2(1,1:ep_idx);
% e22_ep1 = e2(2,1:ep_idx);
% e21_ep2 = e2(1,ep_idx+1:end);
% e22_ep2 = e2(2,ep_idx+1:end);
% e31_ep1 = e3(1,1:ep_idx);
% e32_ep1 = e3(2,1:ep_idx);
% e31_ep2 = e3(1,ep_idx+1:end);
% e32_ep2 = e3(2,ep_idx+1:end);
% e41_ep1 = e4(1,1:ep_idx);
% e42_ep1 = e4(2,1:ep_idx);
% e41_ep2 = e4(1,ep_idx+1:end);
% e42_ep2 = e4(2,ep_idx+1:end);

% ISE = @(e) sqrt(sum(e.^2)*sim_dt);

% fprintf("Norm of error in Episode 1: \n")
% fprintf("C1 e1 ep1: %.3f\n", ISE(e11_ep1)*1e3)
% fprintf("C1 e2 ep1: %.3f\n", ISE(e12_ep1)*1e3)
% fprintf("C2 e1 ep1: %.3f\n", ISE(e21_ep1)*1e3)
% fprintf("C2 e2 ep1: %.3f\n", ISE(e22_ep1)*1e3)
% fprintf("C3 e1 ep1: %.3f\n", ISE(e31_ep1)*1e3)
% fprintf("C3 e2 ep1: %.3f\n", ISE(e32_ep1)*1e3)
% fprintf("C4 e1 ep1: %.3f\n", ISE(e41_ep1)*1e3)
% fprintf("C4 e2 ep1: %.3f\n", ISE(e42_ep1)*1e3)

% fprintf("Norm of error in Episode 2: \n")
% fprintf("C1 e1 ep2: %.3f\n", ISE(e11_ep2)*1e3)
% fprintf("C1 e2 ep2: %.3f\n", ISE(e12_ep2)*1e3)
% fprintf("C2 e1 ep2: %.3f\n", ISE(e21_ep2)*1e3)
% fprintf("C2 e2 ep2: %.3f\n", ISE(e22_ep2)*1e3)
% fprintf("C3 e1 ep2: %.3f\n", ISE(e31_ep2)*1e3)
% fprintf("C3 e2 ep2: %.3f\n", ISE(e32_ep2)*1e3)
% fprintf("C4 e1 ep2: %.3f\n", ISE(e41_ep2)*1e3)
% fprintf("C4 e2 ep2: %.3f\n", ISE(e42_ep2)*1e3)

% fprintf("Improvement in Episode 2: \n")
% fprintf("C1 e1: %.3f\n", 1-ISE(e11_ep2)/ISE(e11_ep1))
% fprintf("C1 e2: %.3f\n", 1-ISE(e12_ep2)/ISE(e12_ep1))
% fprintf("C2 e1: %.3f\n", 1-ISE(e21_ep2)/ISE(e21_ep1))
% fprintf("C2 e2: %.3f\n", 1-ISE(e22_ep2)/ISE(e22_ep1))
% fprintf("C3 e1: %.3f\n", 1-ISE(e31_ep2)/ISE(e31_ep1))
% fprintf("C3 e2: %.3f\n", 1-ISE(e32_ep2)/ISE(e32_ep1))
% fprintf("C4 e1: %.3f\n", 1-ISE(e41_ep2)/ISE(e41_ep1))
% fprintf("C4 e2: %.3f\n", 1-ISE(e42_ep2)/ISE(e42_ep1))

beep()

%% LOCAL FUNCTIONS
function y = time_norm(x)
    num_x = size(x, 1);

    y = zeros(1,size(x,2));
    for x_idx = 1:1:num_x
        y = y + x(x_idx,:).^2;
    end
    y = sqrt(y);
end