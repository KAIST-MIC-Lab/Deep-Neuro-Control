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

% ============================================
%     Fig. 1: State 1 (Ref vs Obs)
% ============================================
figure(1); clf; 
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(t, x_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(t, r_hist(1,:), "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(r_hist(1,:)); minVal = min(r_hist(1,:)); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     Fig. 2: State 2 (Ref vs Obs)
% ============================================
figure(2); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(t, x_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(t, r_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(r_hist(2,:)); minVal = min(r_hist(2,:)); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 3: Control Input 1
% ============================================
figure(3);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];

plot(t, u_hist(1,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(u_hist(1,:)); minVal = min(u_hist(1,:)); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 4: Control Input 2
% ============================================
figure(4);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];

plot(t, u_hist(2,:), "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(u_hist(2,:)); minVal = min(u_hist(2,:)); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';
