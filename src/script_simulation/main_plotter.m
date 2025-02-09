clear

%%
SAVE_FLAG = 1;
POSITION_FLAG = 1; % it will plot fiugures in the same position
gray = "#808080";

more_blue = "#0072BD";
more_red = "#A2142F";

%% 
ctrl1_name = "6-Feb-2025_15-44-57"; % with constraint
ctrl2_name = "6-Feb-2025_15-45-27"; % without constraint

%%
ctrl1_log = load("sim_result/"+ctrl1_name+".mat");
ctrl2_log = load("sim_result/"+ctrl2_name+".mat");

%% RESULT PLOTTER

assert(ctrl1_log.T == ctrl2_log.T, "Simulation time is not the same")
assert(ctrl1_log.u_max == ctrl2_log.u_max, "Control input limit is not the same")

T = ctrl1_log.T;
u_max = ctrl1_log.u_max;
V_max = ctrl1_log.V_max;

ctrl1.id = signal2data(ctrl1_log.logsout, "id");
ctrl1.iq = signal2data(ctrl1_log.logsout, "iq");
ctrl1.id_ref = signal2data(ctrl1_log.logsout, "id_ref");
ctrl1.iq_ref = signal2data(ctrl1_log.logsout, "iq_ref");
ctrl1.vd = signal2data(ctrl1_log.logsout, "vd");
ctrl1.vq = signal2data(ctrl1_log.logsout, "vq");

ctrl1.Vn0 = signal2data(ctrl1_log.logsout, "Vn0");
ctrl1.Vn1 = signal2data(ctrl1_log.logsout, "Vn1");
ctrl1.lbd_th0 = signal2data(ctrl1_log.logsout, "lbd_th0");
ctrl1.lbd_th1 = signal2data(ctrl1_log.logsout, "lbd_th1");
ctrl1.lbd_u = signal2data(ctrl1_log.logsout, "lbd_u");

ctrl2.id = signal2data(ctrl2_log.logsout, "id");
ctrl2.iq = signal2data(ctrl2_log.logsout, "iq");
ctrl2.id_ref = signal2data(ctrl2_log.logsout, "id_ref");
ctrl2.iq_ref = signal2data(ctrl2_log.logsout, "iq_ref");
ctrl2.vd = signal2data(ctrl2_log.logsout, "vd");
ctrl2.vq = signal2data(ctrl2_log.logsout, "vq");

ctrl2.Vn0 = signal2data(ctrl2_log.logsout, "Vn0");
ctrl2.Vn1 = signal2data(ctrl2_log.logsout, "Vn1");
ctrl2.lbd_th0 = signal2data(ctrl2_log.logsout, "lbd_th0");
ctrl2.lbd_th1 = signal2data(ctrl2_log.logsout, "lbd_th1");
ctrl2.lbd_u = signal2data(ctrl2_log.logsout, "lbd_u");

%%
font_size = 24;
line_width = 2;
lgd_size = 16;
    
fig_height = 250; fig_width = 800;

% ============================================
%     Fig. 1: Current d-axis (Ref vs Obs)
% ============================================
figure(1); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

id1 = ctrl1.id; id2 = ctrl2.id; id_ref = ctrl1.id_ref;

plot([0.5 0.5], [-5 5], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(.02, -3.5, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(.52, -3.5, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

plot(id2.Time, id2.Data, "Color", gray, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(id1.Time, id1.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(id_ref.Time, id_ref.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$i_s^d$ / A', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(id_ref.Data); minVal = min(id_ref.Data); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% ============================================
%     Fig. 2: Current q-axis (Ref vs Obs)
% ============================================
figure(2); clf;
hF = gcf; 
hF.Position(3:4) = [fig_width, fig_height];

iq1 = ctrl1.iq; iq2 = ctrl2.iq; iq_ref = ctrl1.iq_ref;

plot([0.5 0.5], [-5 5], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(.02, -3.5, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(.52, -3.5, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

plot(iq2.Time, iq2.Data, "Color", gray, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(iq1.Time, iq1.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(iq_ref.Time, iq_ref.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$i_s^q$ / A', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(iq_ref.Data); minVal = min(iq_ref.Data); 
% maxVal = u_max; minVal = -u_max;
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 3: Voltage Norm
% ============================================
figure(3);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];

plot([0.5 0.5], [-50 50], "Color", "black", "LineWidth", line_width, "LineStyle", "-."); hold on
text(.02, 45, "Episode 1", "FontSize", font_size, "FontName", 'Times New Roman')
text(.52, 45, "Episode 2", "FontSize", font_size, "FontName", 'Times New Roman')

vd1 = ctrl1.vd; vq1 = ctrl1.vq;
norm_v1 = sqrt( vd1.Data.^2 + vq1.Data.^2 );

vd2 = ctrl2.vd; vq2 = ctrl2.vq;
norm_v2 = sqrt( vd2.Data.^2 + vq2.Data.^2 );

plot(vd1.Time, ones(size(vd1.Time))*u_max, "Color", "black", "LineWidth", line_width, "LineStyle", "-"); hold on
plot(vd2.Time, norm_v2, "Color", gray, "LineWidth", line_width, "LineStyle", "-"); hold on
plot(vd1.Time, norm_v1, "Color", "blue", "LineWidth", line_width, "LineStyle", "-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\Vert \mathbf{u}_s^{dq}\Vert$ / V', 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = max(norm_v2); minVal = min(norm_v2); 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 4: Weight
% ============================================
figure(4);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];

Vn0_1 = ctrl1.Vn0; Vn1_1 = ctrl1.Vn1;
Vn0_2 = ctrl2.Vn0; Vn1_2 = ctrl2.Vn1;

plot(Vn0_1.Time, ones(size(Vn0_1.Time))*V_max(1), "Color", "red", "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\bar \theta_0$"); hold on
plot(Vn1_1.Time, ones(size(Vn1_1.Time))*V_max(2), "Color", "blue", "LineWidth", line_width, "LineStyle", "--", "DisplayName", "$\bar \theta_1$"); hold on

plot(Vn0_1.Time, Vn0_1.Data, "Color", more_red, "LineWidth", line_width, "LineStyle", "-", "DisplayName", "With $c_u$: $\Vert \hat\theta_0\Vert$"); hold on
plot(Vn1_1.Time, Vn1_1.Data, "Color", more_blue, "LineWidth", line_width, "LineStyle", "-", "DisplayName", "With $c_u$: $\Vert \hat\theta_1\Vert$"); hold on
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\Vert \theta_i\Vert$', 'FontSize', font_size, 'Interpreter', 'latex');

plot(Vn0_2.Time, Vn0_2.Data, "Color", "red", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "Without $c_u$: $\Vert \hat\theta_1\Vert$"); hold on
plot(Vn1_2.Time, Vn1_2.Data, "Color", "blue", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "Without $c_u$: $\Vert \hat\theta_1\Vert$"); hold on
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    lgd.NumColumns = 3;
    lgd.Location = 'northwest';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 
grid on; grid minor;
maxVal = V_max(2); minVal = 0; 
len = maxVal-minVal; ratio = .1;
ylim([minVal-len*ratio maxVal+len*ratio]);
xlim([0 T])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 5: Multipliers
% ============================================
figure(5);clf
hF = gcf;
hF.Position(3:4) = [fig_width, fig_height];

lbd_th0_1 = ctrl1.lbd_th0; lbd_th1_1 = ctrl1.lbd_th1; lbd_u_1 = ctrl1.lbd_u;
lbd_th0_2 = ctrl2.lbd_th0; lbd_th1_2 = ctrl2.lbd_th1; lbd_u_2 = ctrl2.lbd_u;

semilogy(lbd_th0_1.Time, exp(lbd_th0_1.Data), "Color", more_red, "LineWidth", line_width, "LineStyle", "-", "DisplayName", "With $c_u$: $\lambda_{\theta_0}$"); hold on
semilogy(lbd_th1_1.Time, exp(lbd_th1_1.Data), "Color", more_blue, "LineWidth", line_width, "LineStyle", "-", "DisplayName", "With $c_u$: $\lambda_{\theta_1}$"); hold on

semilogy(lbd_u_1.Time, exp(lbd_u_1.Data), "Color", "green", "LineWidth", line_width, "LineStyle", "-", "DisplayName", "With $c_u$: $\lambda_u$"); hold on

semilogy(lbd_th0_1.Time, exp(lbd_th0_1.Data), "Color", "red", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "Without $c_u$: $\lambda_{\theta_0}$"); hold on
semilogy(lbd_th1_1.Time, exp(lbd_th1_1.Data), "Color", "blue", "LineWidth", line_width, "LineStyle", "-.", "DisplayName", "Without $c_u$: $\lambda_{\theta_1}$"); hold on

xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\lambda_j$ (Log scale)', 'FontSize', font_size, 'Interpreter', 'latex');
    lgd = legend;
    % lgd.Orientation = 'Vertical';
    lgd.NumColumns = 1;
    lgd.Location = 'northwest';
    lgd.Interpreter = 'latex';
    lgd.FontSize = lgd_size; 
grid on; grid minor;
xlim([0 T])
    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';

% ============================================
%        Fig. 6: Model Parameters
% ============================================
plotter_machine_data;

%% SAVE FIGURES
if SAVE_FLAG
    [~,~] = mkdir("figures/compare");

    for idx = 1:1:6

        f_name = "figures/compare/Fig" + string(idx);

        saveas(figure(idx), f_name + ".png")
        exportgraphics(figure(idx), f_name+'.eps')
    end
end

beep()