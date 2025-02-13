% ********************************************
%
% Change the name of the controller to the one you want to plot. (ctrl_name)
%
% ********************************************

%% FASTEN YOUR SEATBELT
clear

%% FIGURE SETTING
RESULT_SAVE_FLAG = 1;
POSITION_FLAG = 1; % it will plot fiugures in the same position
FIGURE_SAVE_FLAG = 1;

font_size = 16;
line_width = 2;
lgd_size = 16;
% font_size = 32;
% line_width = 2;
% lgd_size = 28;

fig_height = 200; 
fig_width = 450;
% For papers
% fig_height = 300; 
% fig_width = 800;

%% DATA LOAD
ctrl_name = "9-Feb-2025_23-25-01";
ctrl_hist = load("results/"+ctrl_name+".mat");

%% DATA PROCESSING
T = ctrl_hist.T;

x1 = signal2data(ctrl_hist.logsout, "x1");
x2 = signal2data(ctrl_hist.logsout, "x2");
u1 = signal2data(ctrl_hist.logsout, "u1");
u2 = signal2data(ctrl_hist.logsout, "u2");
r1 = signal2data(ctrl_hist.logsout, "r1");
r2 = signal2data(ctrl_hist.logsout, "r2");

%% MAIN PLOT FUNCTIONS

% ============================================
%     Fig. 1: State 1 (Ref vs Obs)
% ============================================
figure(1); clf; 
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(r1.Time, r1.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(x1.Time, x1.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(r1.Data); minVal = min(r1.Data); 
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

plot(r2.Time, r2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(x2.Time, x2.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(r2.Data); minVal = min(r2.Data); 
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

plot(u1.Time, u1.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(u1.Data); minVal = min(u1.Data); 
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

plot(u2.Time, u2.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(u2.Data); minVal = min(u2.Data); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';



%% FIGURE SAVE
if FIGURE_SAVE_FLAG
    [~,~] = mkdir("figures/"+ctrl_name);

    for idx = 1:1:4   
        f_name = "figures/" + ctrl_name + "/Fig" + string(idx);

        saveas(figure(idx), f_name + ".png")
        exportgraphics(figure(idx), f_name+'.eps')
    end
end
