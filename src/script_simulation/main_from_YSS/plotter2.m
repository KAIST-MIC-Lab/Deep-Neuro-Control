%% FIGURE SETTING
POSITION_FLAG = 1; % it will plot figures in the same position

font_size = 16;
line_width = 2;
lgd_size = 16;
fig_height = 200;
fig_width  = 450;

% For papers
% font_size = 32;
% line_width = 2;
% lgd_size = 28;
% fig_height = 300;
% fig_width  = 800;

%% MAIN PLOT FUNCTIONS
figure(1); clf;
hF = gcf;
% 기본 크기
hF.Position(3:4) = [1600, 800];

% 창을 화면에서 '조금 아래'로 내리기 (모니터에서 위쪽에 뜨는 문제 완화)
try
    hF.Position(2) = max(0, hF.Position(2) - 120); % 아래로 120px
catch
    % Position 설정 불가한 환경이면 무시
end

tl = tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     Fig. 1: State(Ref vs Obs) - x1
% ============================================
nexttile;
plot(t, x_hist(1,:),    "Color","blue",   "LineWidth",line_width, "LineStyle","-"); hold on
plot(t, x_non_hist(1,:), "Color","cyan",   "LineWidth",line_width, "LineStyle","-");
plot(t, xd_hist(1,:),   "Color","red",    "LineWidth",line_width, "LineStyle","--");
plot(t, th_ref_hist(1,:), "Color","magenta","LineWidth",line_width, "LineStyle","--");

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$x_1$',    'FontSize', font_size, 'Interpreter','latex');

maxVal = max(xd_hist(1,:)); minVal = min(xd_hist(1,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$x_1$ (with disturbance)', ...
        '$x_1$ (no disturbance)', ...
        '$x_{1,\mathrm{des}}$', ...
        '$\theta_{\mathrm{ref}}$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',2);

% ============================================
%     Fig. 1: State(Ref vs Obs) - x2
% ============================================
nexttile;
plot(t, x_hist(2,:),    "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on
plot(t, x_non_hist(2,:), "Color","cyan", "LineWidth",line_width, "LineStyle","-");
plot(t, xd_hist(2,:),   "Color","red",  "LineWidth",line_width, "LineStyle","--");

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$x_2$',    'FontSize', font_size, 'Interpreter','latex');

maxVal = max(xd_hist(2,:)); minVal = min(xd_hist(2,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$x_2$ (with disturbance)', ...
        '$x_2$ (no disturbance)', ...
        '$x_{2,\mathrm{des}}$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',2);

% ============================================
%     Fig. 2: Condition number (or metric norm) M_hist
% ============================================
nexttile;
plot(t, M_hist(1,:), "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$\|M^{-1}\|$', 'FontSize', font_size, 'Interpreter','latex');

maxVal = max(M_hist(1,:)); minVal = min(M_hist(1,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$\| (W_{\mathrm{bar}}/\mu)^{-1} \|$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',1);

% ============================================
%     Fig. 2: Optimization done flag
% ============================================
nexttile;
plot(t, optDone_hist(1,:), "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('\texttt{optDone}', 'FontSize', font_size, 'Interpreter','latex');
xlim([0 T]);
ylim([-0.1 1.1]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'done flag'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',1);

% ============================================
%     Fig. 2: mu history
% ============================================
nexttile;
plot(t, mu_hist(1,:), "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$\mu$',       'FontSize', font_size, 'Interpreter','latex');

maxVal = max(mu_hist(1,:)); minVal = min(mu_hist(1,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$\mu$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',1);

% ============================================
%     Fig. 2: Control input u1 vs ud1
% ============================================
nexttile;
plot(t, u_hist(1,:), "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on
plot(t, ud_hist(1,:), "Color","red",  "LineWidth",line_width, "LineStyle","--");

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$u_1$',      'FontSize', font_size, 'Interpreter','latex');

maxVal = max(ud_hist(1,:)); minVal = min(ud_hist(1,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$u_1$', '$u_{1,\mathrm{des}}$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',2);

% ============================================
%     Fig. 2: Control input u2 vs ud2
% ============================================
nexttile;
plot(t, u_hist(2,:), "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on
plot(t, ud_hist(2,:), "Color","red",  "LineWidth",line_width, "LineStyle","--");

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$u_2$',      'FontSize', font_size, 'Interpreter','latex');

maxVal = max(ud_hist(2,:)); minVal = min(ud_hist(2,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$u_2$', '$u_{2,\mathrm{des}}$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',2);

% ============================================
%     Fig. 2: Control input norm
% ============================================
nexttile;
u_norm = sqrt(sum(u_hist.^2));

plot(t, u_norm, "Color","blue", "LineWidth",line_width, "LineStyle","-"); hold on
% plot(t, ud_hist(2,:), "Color","red",  "LineWidth",line_width, "LineStyle","--");

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter','latex');
ylabel('$\Vert u\Vert$',      'FontSize', font_size, 'Interpreter','latex');

maxVal = max(u_norm(2,:)); minVal = min(u_norm(2,:));
len = maxVal - minVal; ratio = .1;
if len ~= 0
    ylim([minVal - len*ratio, maxVal + len*ratio]);
end
xlim([0 T]);

ax = gca; ax.FontSize = font_size; ax.FontName = 'Times New Roman';

legend({'$u_2$', '$u_{2,\mathrm{des}}$'}, ...
       'Interpreter','latex','FontSize',lgd_size, ...
       'Location','southoutside','NumColumns',2);

% ============================================
% (빈 타일 정리: 필요시 사용)
% ============================================
% nexttile; axis off
% nexttile; axis off

% 전체 타이틀이 필요하면 사용
% title(tl, 'NCM-based Control Results', 'Interpreter','latex', 'FontSize', font_size+2);
