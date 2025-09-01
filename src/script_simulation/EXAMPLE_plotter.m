clear

%%
SAVE_FLAG = 1;
POSITION_FLAG = 1; % it will plot fiugures in the same position
TEX_CONVERT_FLAG = 0;
gray = "#808080";

more_blue = "#0072BD";
more_red = "#A2142F";


%%
if TEX_CONVERT_FLAG
    font_size = 4;
    line_width = .5;
    lgd_size = 2;

    fig_height = 8.89/4;
    fig_width = 8.89;

    fig_unit = 'centimeters';
else
    font_size = 17;
    % font_size = 17;
    line_width = 1.5;
    lgd_size = 16;
        
    fig_height = 200;
    % fig_height = 150;
    fig_width = 800;

    fig_unit = 'pixels';
end

%% 
% ctrl1_name = "c1.trc"; 
ctrl2_name = "c2.trc";
ctrl3_name = "c3.trc"; 
ctrl4_name = "c4.trc"; 

ctrl1_name = "re_c1.trc"; 
% ctrl2_name = "re_c3.trc";
% ctrl3_name = "re_c4.trc"; 
% ctrl4_name = "re_c2.trc"; 

data1 = trc2data(ctrl1_name, 1);
data2 = trc2data(ctrl2_name, 2);
data3 = trc2data(ctrl3_name, 3);
data4 = trc2data(ctrl4_name, 4);


c1 = "blue"; c2="cyan"; c3=gray; c4="magenta";
th_max = [6 6 6];
u_ball = 11;
u_max2 = 3.5;
u_max1 = sqrt(u_ball^2 - u_max2^2);


%%
idleTime1 = 3+8;
idleTime2 = 4;

ep_time = 12; T = ep_time * 2;

start_t = 29.4;
end_t   = 31.9;

% general_lim = [start_t end_t];
general_lim = [idleTime1 idleTime1+T];

focus_c1 = find(data1.u1.Time >= start_t & data1.u1.Time <= end_t);
focus_c2 = find(data2.u1.Time >= start_t & data2.u1.Time <= end_t);
focus_c3 = find(data3.u1.Time >= start_t & data3.u1.Time <= end_t);
focus_c4 = find(data4.u1.Time >= start_t & data4.u1.Time <= end_t);

%% FIG. 1; STATE 1
figure(1); clf;
hF = gcf; 
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.r1.Time, data4.r1.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data3.r1.Time, data3.r1.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data2.r1.Time, data2.r1.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data1.r1.Time, data1.r1.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

plot(data4.q1.Time, data4.q1.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.q1.Time, data3.q1.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.q1.Time, data2.q1.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.q1.Time, data1.q1.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(idleTime1+.2,         -1.3, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(idleTime1+ep_time+.2, -1.3, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$q_1$ / rad', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(data1.r1.Data); minVal = min(data1.r1.Data); 
maxVal = 1.5; minVal = -1.5;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 2; STATE 2
figure(2); clf;
hF = gcf; 
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.r2.Time, data4.r2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data3.r2.Time, data3.r2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data2.r2.Time, data2.r2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
plot(data1.r2.Time, data1.r2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

plot(data4.q2.Time, data4.q2.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.q2.Time, data3.q2.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.q2.Time, data2.q2.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.q2.Time, data1.q2.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(idleTime1+.2,         -2.3, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(idleTime1+ep_time+.2, -2.3, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$q_2$ / rad', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = max(data1.r2.Data); minVal = min(data1.r2.Data); 
maxVal = 1; minVal = -2.5;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 3; CONTROL INPUT 1
figure(3);clf
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.u1.Time, data4.u1.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.u1.Time, data3.u1.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.u1.Time, data2.u1.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.u1.Time, data1.u1.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

plot([-10 T+15], +1*[u_max1 u_max1], "Color", "black", "LineWidth", line_width, "LineStyle", "--"); hold on
plot([-10 T+15], -1*[u_max1 u_max1], "Color", "black", "LineWidth", line_width, "LineStyle", "--"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(idleTime1+.2,         2.3, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(idleTime1+ep_time+.2, 2.3, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\tau_1$ / Nm', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = 14; minVal = 1;
% maxVal = max(norm_v2); minVal = min(norm_v2); 
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 4; CONTROL INPUT 2
figure(4);clf
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.u2.Time, data4.u2.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.u2.Time, data3.u2.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.u2.Time, data2.u2.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.u2.Time, data1.u2.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

plot([-10 T+15], +1*[u_max2 u_max2], "Color", "black", "LineWidth", line_width, "LineStyle", "--"); hold on
plot([-10 T+15], -1*[u_max2 u_max2], "Color", "black", "LineWidth", line_width, "LineStyle", "--"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(idleTime1+.2,         -.65, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(idleTime1+ep_time+.2, -.65, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\tau_2$ / Nm', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = 4; minVal = -1;
% maxVal = max(norm_v2); minVal = min(norm_v2); 
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 5; CONTROL INPUT NORM
figure(5);clf
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

norm_u1 = sqrt(data1.u1.Data.^2 + data1.u2.Data.^2);
norm_u2 = sqrt(data2.u1.Data.^2 + data2.u2.Data.^2);
norm_u3 = sqrt(data3.u1.Data.^2 + data3.u2.Data.^2);
norm_u4 = sqrt(data4.u1.Data.^2 + data4.u2.Data.^2);

plot(data4.u1.Time, norm_u4, "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.u1.Time, norm_u3, "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.u1.Time, norm_u2, "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.u1.Time, norm_u1, "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

plot([-10 T+15], [u_ball u_ball], "Color", "black", "LineWidth", line_width, "LineStyle", "--"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(idleTime1+.2,         2.5, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(idleTime1+ep_time+.2, 2.5, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\Vert {\mbox{\boldmath $\tau$}} \Vert$ / Nm', 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = u_ball; minVal = 0; 
maxVal = 13; minVal = 2;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 6; CONTROL NORM (BIRDEYE)
figure(6);clf
hf = gcf;
hF.Units = fig_unit;
hf.Position(3:4) = [fig_width, fig_height];

ang = 0:0.01:2*pi;

plot(u_ball*cos(ang), u_ball*sin(ang), "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-100, 100], [1, 1] * u_max2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-100, 100], [-1, -1] * u_max2, "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([1, 1] * u_max1, [-100, 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on
plot([-1, -1] * u_max1, [-100, 100], "color", 'black', "LineWidth", line_width, "LineStyle", "-."); hold on

plot(data4.u1.Data(focus_c1), data4.u2.Data(focus_c1), "Color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.u1.Data(focus_c2), data3.u2.Data(focus_c2), "Color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.u1.Data(focus_c3), data2.u2.Data(focus_c3), "Color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.u1.Data(focus_c4), data1.u2.Data(focus_c4), "Color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

xlabel("$\tau_1$ / Nm", "Interpreter", "latex")
ylabel("$\tau_2$ / Nm", "Interpreter", "latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on; grid minor;

xlim([10 13])
ylim([1.5 3.7])

% pbaspect([1 1 1])

%% FIG. 7; WEIGHT NORM
figure(7);clf
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

% plot([-10 T+100], +1*[th_max(1) th_max(1)], "Color", "black", "LineWidth", line_width, "LineStyle", "--", 'HandleVisibility','off'); hold on
% plot([-10 T+100], +1*[th_max(2) th_max(2)], "Color", "black", "LineWidth", line_width, "LineStyle", "-.", 'HandleVisibility','off'); hold on
plot([-10 T+100], +1*[th_max(3) th_max(3)], "Color", "black", "LineWidth", line_width, "LineStyle", "--",  'HandleVisibility','off'); hold on

plot(data1.th0.Time, data1.th0.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
plot(data1.th1.Time, data1.th1.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data1.th2.Time, data1.th2.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

plot(data2.th0.Time, data2.th0.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
plot(data2.th1.Time, data2.th1.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data2.th2.Time, data2.th2.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

plot(data3.th0.Time, data3.th0.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
plot(data3.th1.Time, data3.th1.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data3.th2.Time, data3.th2.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

plot(data4.th0.Time, data4.th0.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
plot(data4.th1.Time, data4.th1.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data4.th2.Time, data4.th2.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on
text(idleTime1+.2,         1.85, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman', "HandleVisibility", "off") 
text(idleTime1+ep_time+.2, 1.85, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman', "HandleVisibility", "off")

xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\Vert \widehat{\mbox{\boldmath $\theta$}}_i \Vert$', 'FontSize', font_size, 'Interpreter', 'latex');
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    % lgd.Orientation = 'Horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    % lgd.FontSize = lgd_size; 
    lgd.FontSize = 12; 
    
    
grid on; grid minor;
maxVal = max(th_max); minVal = 1.2; 
len = maxVal-minVal; ratio = .1;
ylim([minVal maxVal+len*ratio]);
xlim(general_lim)
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 8; LAGRANGE MULTIPLIERS
figure(8);clf
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

semilogy(    data1.lbdu.Time,     data1.lbdu.Data, "Color", c1, "LineWidth", line_width, "LineStyle",  "-", "DisplayName", "$\lambda_{{\mbox{\boldmath $\tau$}}}$"); hold on
semilogy(data1.lbdu2Max.Time, data1.lbdu2Max.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\lambda_{\overline{\tau}_2}$"); hold on
semilogy(data1.lbdu2Min.Time, data1.lbdu2Min.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\lambda_{\underline{\tau}_2}$"); hold on

semilogy(    data2.lbdu.Time,     data2.lbdu.Data, "Color", c2, "LineWidth", line_width, "LineStyle",  "-", "DisplayName", "$\lambda_{{\mbox{\boldmath $\tau$}}}$"); hold on
semilogy(data2.lbdu2Max.Time, data2.lbdu2Max.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\lambda_{\overline{\tau}_2}$"); hold on
semilogy(data2.lbdu2Min.Time, data2.lbdu2Min.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\lambda_{\underline{\tau}_2}$"); hold on

semilogy(data3.lbdth2.Time,   data3.lbdth2.Data,   "Color", c3, "LineWidth", line_width, "LineStyle", "-", "DisplayName", "$\lambda_{\theta_2}$"); hold on

xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\lambda_j$ (Log scale)', 'FontSize', font_size, 'Interpreter', 'latex');
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    lgd.NumColumns = 3;
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

grid on; grid minor;
% maxVal = max(data1.lbdu.Data); minVal = 0;
maxVal = 5e-1; minVal = 1e-4;
% len = maxVal-minVal; ratio = .3;
% ylim([minVal-len*ratio maxVal+len*ratio]);
ylim([minVal maxVal]);
xlim([start_t end_t])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

%% FIG. 9; AUXILIARY STATE
figure(9); clf;
hF = gcf; 
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.z1.Time, data4.z1.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-", 'DisplayName', '$\zeta_1$'); hold on
plot(data4.z2.Time, data4.z2.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-.",'DisplayName', '$\zeta_2$'); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\zeta_1,\zeta_2$', 'FontSize', font_size, 'Interpreter', 'latex');
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    % lgd.NumColumns = 3;
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 

% maxVal = max(data2.zeta(1,data2.obs)); minVal = 0; 
maxVal = 3.5; minVal = 0; 
len = maxVal-minVal; ratio = .3;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([start_t end_t])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% =============================================
%    Fig. 10: Computational Time
% =============================================
figure(10); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

plot(data4.cmp.Time, data4.cmp.Data, "color", c4, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data3.cmp.Time, data3.cmp.Data, "color", c3, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data2.cmp.Time, data2.cmp.Data, "color", c2, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(data1.cmp.Time, data1.cmp.Data, "color", c1, "LineWidth", line_width, "LineStyle", "-"); hold on

xlabel("Time / s", "Interpreter", "latex")
% ylabel("Comp. Time / ms", "Interpreter", "latex")
ylabel("$T_{\rm{comp}}$ / ms", "Interpreter", "latex")
set(gca, 'FontSize', font_size, 'FontName', 'Times New Roman')
grid on 
xlim(general_lim)
ylim([3.1 3.7])

%% Fig. 11: (Scope) weight norms
figure(11); clf;
hF = gcf;
hF.Units = fig_unit;
hF.Position(3:4) = [fig_width, fig_height];


% plot([-10 T+100], +1*[th_max(1) th_max(1)], "Color", "black", "LineWidth", line_width, "LineStyle", "--", 'HandleVisibility','off'); hold on
% plot([-10 T+100], +1*[th_max(2) th_max(2)], "Color", "black", "LineWidth", line_width, "LineStyle", "-.", 'HandleVisibility','off'); hold on
plot([-10 T+100], +1*[th_max(3) th_max(3)], "Color", "black", "LineWidth", line_width, "LineStyle", "--",  'HandleVisibility','off'); hold on

% plot(data1.th0.Time, data1.th0.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
% plot(data1.th1.Time, data1.th1.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data1.th2.Time, data1.th2.Data, "Color", c1, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

% plot(data2.th0.Time, data2.th0.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
% plot(data2.th1.Time, data2.th1.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data2.th2.Time, data2.th2.Data, "Color", c2, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

% plot(data3.th0.Time, data3.th0.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
% plot(data3.th1.Time, data3.th1.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data3.th2.Time, data3.th2.Data, "Color", c3, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

% plot(data4.th0.Time, data4.th0.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_0$"); hold on
% plot(data4.th1.Time, data4.th1.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_1$"); hold on
plot(data4.th2.Time, data4.th2.Data, "Color", c4, "LineWidth", line_width, "LineStyle", "-",  "DisplayName", "$\widehat{\mbox{\boldmath $\theta$}}_2$"); hold on

plot(idleTime1+[ep_time ep_time], [-5e1 5e1], "Color", "black", "LineWidth", line_width, "LineStyle", "-.", "HandleVisibility", "off"); hold on
text(idleTime1+.2,         1.6, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman', "HandleVisibility", "off") 
text(idleTime1+ep_time+.2, 1.6, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman', "HandleVisibility", "off")

xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\Vert \widehat{\mbox{\boldmath $\theta$}}_i \Vert$', 'FontSize', font_size, 'Interpreter', 'latex');
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    % lgd.Orientation = 'Horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southeast';
    lgd.Interpreter = 'latex';
    % lgd.FontSize = lgd_size; 
    lgd.FontSize = 12; 
    
    
grid on; grid minor;
% maxVal = max(th_max); minVal = 1.2; 
maxVal = 6.2; minVal = 5.2;
len = maxVal-minVal; ratio = .1;
ylim([minVal maxVal+len*ratio]);
xlim([start_t end_t])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';


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
ctrl_dt = 1/250;
sim_dt = ctrl_dt / 1000;

e1 = [
    transpose(data1.q1.Data-data1.r1.Data);
    transpose(data1.q2.Data-data1.r2.Data)
];
e1 = e1(:, (data1.r1.Time >= idleTime1 & data1.r1.Time <= idleTime1+T));
e2 = [
    transpose(data2.q1.Data-data2.r1.Data);
    transpose(data2.q2.Data-data2.r2.Data)
];
e2 = e2(:, (data2.r1.Time >= idleTime1 & data2.r1.Time <= idleTime1+T));
e3 = [
    transpose(data3.q1.Data-data3.r1.Data);
    transpose(data3.q2.Data-data3.r2.Data)
];
e3 = e3(:, (data3.r1.Time >= idleTime1 & data3.r1.Time <= idleTime1+T));
e4 = [
    transpose(data4.q1.Data-data4.r1.Data);
    transpose(data4.q2.Data-data4.r2.Data)
];
e4 = e4(:, (data4.r1.Time >= idleTime1 & data4.r1.Time <= idleTime1+T));

ep_idx = floor(size(e1,2)/2);

e11_ep1 = e1(1,1:ep_idx);
e12_ep1 = e1(2,1:ep_idx);
e11_ep2 = e1(1,ep_idx+1:end);
e12_ep2 = e1(2,ep_idx+1:end);
e21_ep1 = e2(1,1:ep_idx);
e22_ep1 = e2(2,1:ep_idx);
e21_ep2 = e2(1,ep_idx+1:end);
e22_ep2 = e2(2,ep_idx+1:end);
e31_ep1 = e3(1,1:ep_idx);
e32_ep1 = e3(2,1:ep_idx);
e31_ep2 = e3(1,ep_idx+1:end);
e32_ep2 = e3(2,ep_idx+1:end);
e41_ep1 = e4(1,1:ep_idx);
e42_ep1 = e4(2,1:ep_idx);
e41_ep2 = e4(1,ep_idx+1:end);
e42_ep2 = e4(2,ep_idx+1:end);

ISE = @(e) sqrt(sum(e.^2)*sim_dt);

fprintf("Norm of error in Episode 1: \n")
fprintf("C1 e1 ep1: %.3f\n", ISE(e11_ep1)*1e3)
fprintf("C1 e2 ep1: %.3f\n", ISE(e12_ep1)*1e3)
fprintf("C2 e1 ep1: %.3f\n", ISE(e21_ep1)*1e3)
fprintf("C2 e2 ep1: %.3f\n", ISE(e22_ep1)*1e3)
fprintf("C3 e1 ep1: %.3f\n", ISE(e31_ep1)*1e3)
fprintf("C3 e2 ep1: %.3f\n", ISE(e32_ep1)*1e3)
fprintf("C4 e1 ep1: %.3f\n", ISE(e41_ep1)*1e3)
fprintf("C4 e2 ep1: %.3f\n", ISE(e42_ep1)*1e3)

fprintf("Norm of error in Episode 2: \n")
fprintf("C1 e1 ep2: %.3f\n", ISE(e11_ep2)*1e3)
fprintf("C1 e2 ep2: %.3f\n", ISE(e12_ep2)*1e3)
fprintf("C2 e1 ep2: %.3f\n", ISE(e21_ep2)*1e3)
fprintf("C2 e2 ep2: %.3f\n", ISE(e22_ep2)*1e3)
fprintf("C3 e1 ep2: %.3f\n", ISE(e31_ep2)*1e3)
fprintf("C3 e2 ep2: %.3f\n", ISE(e32_ep2)*1e3)
fprintf("C4 e1 ep2: %.3f\n", ISE(e41_ep2)*1e3)
fprintf("C4 e2 ep2: %.3f\n", ISE(e42_ep2)*1e3)

fprintf("Improvement in Episode 2: \n")
fprintf("C1 e1: %.3f\n", 1-ISE(e11_ep2)/ISE(e11_ep1))
fprintf("C1 e2: %.3f\n", 1-ISE(e12_ep2)/ISE(e12_ep1))
fprintf("C2 e1: %.3f\n", 1-ISE(e21_ep2)/ISE(e21_ep1))
fprintf("C2 e2: %.3f\n", 1-ISE(e22_ep2)/ISE(e22_ep1))
fprintf("C3 e1: %.3f\n", 1-ISE(e31_ep2)/ISE(e31_ep1))
fprintf("C3 e2: %.3f\n", 1-ISE(e32_ep2)/ISE(e32_ep1))
fprintf("C4 e1: %.3f\n", 1-ISE(e41_ep2)/ISE(e41_ep1))
fprintf("C4 e2: %.3f\n", 1-ISE(e42_ep2)/ISE(e42_ep1))

beep()