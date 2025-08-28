% ********************************************
%
% This script is called by main.m to plot the simulation results.
%
% ********************************************

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

%% MAIN PLOT FUNCTIONS

figure(1); clf; 
% hF = gcf; 
% hF.Position(3:4) = [1600, 800];
tl = tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact');


% ============================================
%     State(Ref vs Obs) 1
% ============================================
nexttile;

for idx = 1:1:num_sample
    plot(t, sys{idx}.x_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, sys{idx}.x_non_hist(1,:), "Color", "cyan", "LineWidth", line_width, "LineStyle", "-"); hold on
end
plot(t, xd_hist(1,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(xd_hist(1,:)); minVal = min(xd_hist(1,:)); 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     State(Ref vs Obs) 2
% ============================================
nexttile;
for idx = 1:1:num_sample
    plot(t, sys{idx}.x_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(t, sys{idx}.x_non_hist(2,:), "Color", "cyan", "LineWidth", line_width, "LineStyle", "-"); hold on
end
plot(t, xd_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(xd_hist(2,:)); minVal = min(xd_hist(2,:)); 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     State(Ref vs Obs) 3
% ============================================
% nexttile;
% plot(t, x_hist(3,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
% plot(t, x_non_hist(3,:), "Color", "cyan", "LineWidth", line_width, "LineStyle", "-"); hold on
% plot(t, xd_hist(3,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$x_3$', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(xd_hist(3,:)); minVal = min(xd_hist(3,:)); 
% len = maxVal-minVal; ratio = .1;
% ylim([minVal-len*ratio maxVal+len*ratio]);
% xlim([0 T])

% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';

% ============================================
%     State(Ref vs Obs) 4
% ============================================
% nexttile;
% plot(t, x_hist(4,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
% plot(t, x_non_hist(4,:), "Color", "cyan", "LineWidth", line_width, "LineStyle", "-"); hold on
% plot(t, xd_hist(4,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$x_4$', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(xd_hist(4,:)); minVal = min(xd_hist(4,:)); 
% len = maxVal-minVal; ratio = .1;
% ylim([minVal-len*ratio maxVal+len*ratio]);
% xlim([0 T])

% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';

% ============================================
%     Control input 1
% ============================================
nexttile;
for idx = 1:1:num_sample
    plot(t, sys{idx}.u_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
end
plot(t, ud_hist(1,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(ud_hist(1,:)); minVal = min(ud_hist(1,:)); 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     Control input 2
% ============================================
nexttile;
for idx = 1:1:num_sample
    plot(t, sys{idx}.u_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "--"); hold on
    plot(t, sys{idx}.uSat_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
end
plot(t, ud_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(ud_hist(2,:)); minVal = min(ud_hist(2,:)); 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     Control input loci
% ============================================
nexttile;
for idx = 1:1:num_sample
    plot(sys{idx}.u_hist(1,:), sys{idx}.u_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
end
rad = -pi:0.1:pi;
plot(sin(rad), cos(rad), "Color", "black")

grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(xd_hist(2,:)); minVal = min(xd_hist(2,:)); 
maxVal = max_u; minVal = -max_u;
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
    xlim([minVal-len*ratio maxVal+len*ratio]);
end

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     X
% ============================================
nexttile;
% plot(t, sys{idx}.X_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
for idx = 1:1:num_sample
    plot(t, sys{idx}.M_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
end

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$X$', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(X_hist(1,:)); minVal = min(X_hist(1,:)); 
% maxVal = max(M_hist(1,:)); minVal = min(M_hist(1,:)); 
% len = maxVal-minVal; ratio = .1;
% if len ~= 0
%     ylim([minVal-len*ratio maxVal+len*ratio]);
% end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     optDone
% ============================================
nexttile;
plot(t, sys{idx}.optDone_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Done', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(optDone_hist(1,:)); minVal = min(optDone_hist(1,:)); 
maxVal = 1; minVal = 0; 
len = maxVal-minVal; ratio = .1;
% if len ~= 0
%     ylim([minVal-len*ratio maxVal+len*ratio]);
% end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     mu
% ============================================
nexttile;
plot(t, sys{idx}.mu_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('mu', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(optDone_hist(1,:)); minVal = min(optDone_hist(1,:)); 
maxVal = 1; minVal = 0; 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';


%%
figure(2); clf; 
% hF = gcf; 
% hF.Position(3:4) = [1600, 800];
tl = tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact');

nexttile;
% ============================================
%     Tracking error 1
% ============================================

for idx = 1:1:num_sample
    err = sys{idx}.x_hist(1,:)-xd_hist(1,:);
    plot(t, err, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
end

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(err); minVal = min(err); 
len = maxVal-minVal; ratio = .1;
if len ~= 0
    ylim([minVal-len*ratio maxVal+len*ratio]);
end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';