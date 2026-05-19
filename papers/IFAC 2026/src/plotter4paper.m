

%% FIGURE SETTING
FIGURE_SAVE_FLAG = 1;   % save the figure as .png and .eps

font_size = 22;
line_width = 2.5;
lgd_size = 2;

fig_height = 180;
fig_width = 800/2;

x_list = ["x", "y", "z"];
u_list = ["u_x", "u_y", "u_z"];
x_vec = "\xi";

%% MAIN PLOT

%% MAIN PLOT FUNCTIONS
figure(1); clf; 
hF = gcf; 
hF.Position(3:4) = [800, 600];
% tl = tiledlayout(3, 2, 'Padding', 'none', 'TileSpacing', 'compact');
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'tight');

% f_idx = 1;

% ============================================
%     STATE VARIABLE
% ============================================
for x_idx = 1:1:num_x
    % figure(f_idx); clf; 
    % f_idx = f_idx + 1;
    % hF = gcf; 
    % hF.Position(3:4) = [fig_width, fig_height];
    nexttile(2*(x_idx-1)+1);

    maxVal = -9999; minVal = 9999;
    for c_idx = 1:1:case_num
        plot(t, recs{c_idx}.x_hist(x_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

        maxVal = max(maxVal, max([recs{c_idx}.x_hist(x_idx,:), xd_hist(x_idx,:)]));
        minVal = min(minVal, min([recs{c_idx}.x_hist(x_idx,:), xd_hist(x_idx,:)]));
    end
    plot(t, xd_hist(x_idx,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

    grid on; grid minor;
    xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
    ylabel("$"+x_list(x_idx)+"$", 'FontSize', font_size, 'Interpreter', 'latex');
    len = maxVal-minVal; ratio = .1;
    if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
    xlim([0 T])

    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';
end

% ============================================
%     CONTROL INPUT
% ============================================
for u_idx = 1:1:num_u
    % figure(f_idx); clf; 
    % f_idx = f_idx + 1;
    % hF = gcf; 
    % hF.Position(3:4) = [fig_width, fig_height];
    nexttile(2*(u_idx));

    plot(t, recs{2}.uSat_hist(u_idx,:), "Color", "#008B8B", "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = -9999; minVal = 9999;
    for c_idx = 1:1:case_num
        plot(t, recs{c_idx}.u_hist(u_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
        % plot(t, recs{c_idx}.uSat_hist(u_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "--"); hold on

        maxVal = max(maxVal, max([recs{c_idx}.u_hist(u_idx,:), recs{c_idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
        minVal = min(minVal, min([recs{c_idx}.u_hist(u_idx,:), recs{c_idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
    end

    plot(t, ud_hist(u_idx,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

    grid on; grid minor;
    xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
    ylabel("$"+u_list(u_idx)+"$", 'FontSize', font_size, 'Interpreter', 'latex');
    len = maxVal-minVal; ratio = .1;
    if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
    xlim([0 T])

    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';
end

%% 
figure(2); clf; 
hF = gcf; 
hF.Position(3:4) = [800, 600];
% tl = tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'compact');
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'tight');

% ============================================
%     Control norm
% ============================================
nexttile;
% figure(f_idx); clf; 
% f_idx = f_idx + 1;
% hF = gcf; 
% hF.Position(3:4) = [fig_width, fig_height];

plot([-10, T+10], param.max_u * [1 1], "Color", "black", "LineWidth", line_width, "LineStyle", "-"); hold on
maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    u_norm = time_norm(recs{c_idx}.u_hist);
    plot(t, u_norm, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(u_norm));
    minVal = min(minVal, min(u_norm));
end
text(.3+.02, param.max_u+70, "Input limit $\overline{u}$", "FontSize", font_size-4, "FontName", 'Times New Roman', "Interpreter", "latex")

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel("$\Vert {\mbox{\boldmath $u$}} \Vert$", 'FontSize', font_size, 'Interpreter', 'latex');
% maxVal = param.max_u; minVal = 0;
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     Tracking error 
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num

    err = time_norm(recs{c_idx}.x_hist-xd_hist);
    % err = recs{c_idx}.x_hist(x_idx,:) - xd_hist(x_idx,:);
    plot(t, err, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(err));
    minVal = min(minVal, min(err));
end

grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel("$\Vert {\mbox{\boldmath $e$}}\Vert$", 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     X
% ============================================
nexttile;
% figure(f_idx); clf; 
% f_idx = f_idx + 1;
% hF = gcf; 
% hF.Position(3:4) = [fig_width, fig_height];


maxVal = -9999; minVal = 9999;
warning('X plotted only c_idx \in {2,3} ')
for c_idx = [2,3]
    plot(t, recs{c_idx}.X_hist, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    
    maxVal = max(maxVal, max(recs{c_idx}.X_hist));
    minVal = min(minVal, min(recs{c_idx}.X_hist));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\chi$', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ==============================================
%    nu
% ==============================================
nexttile;

maxVal = -9999; minVal = 9999;
warning('nu plotted only c_idx \in {2} ')
for c_idx = 2
    plot(t, recs{c_idx}.nu_hist, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    
    maxVal = max(maxVal, max(recs{c_idx}.nu_hist));
    minVal = min(minVal, min(recs{c_idx}.nu_hist));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Scaling factor $\nu$', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('Penalty term $\nu$', 'FontSize', font_size, 'Interpreter', 'latex');
minVal = 50;
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%    slack variable (stored in nu)
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
warning('slack variable plotted only c_idx \in {3} ')
for c_idx = 3
    plot(t, recs{c_idx}.nu_hist, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(recs{c_idx}.nu_hist));
    minVal = min(minVal, min(recs{c_idx}.nu_hist));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Slack var. $s$', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     contraction condition satisfaction
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
warning('nu plotted only c_idx \in {2,3} ')
for c_idx = [2,3]
    plot(t, recs{c_idx}.contraction_flag_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    maxVal = 1; minVal = 0;
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('(11) satisfied', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('Contraction Cond.', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

yticks([0 1])
yticklabels({'false','true'})

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';



%% SAVE FIGURES
if FIGURE_SAVE_FLAG
    [~,~] = mkdir("figures/compare");

    for idx = 1:1:2

        f_name = "figures/compare/Fig" + string(idx);

        saveas(figure(idx), f_name + ".png")

        figure(idx);
        % set(gcf, 'Position', [0, 0, fig_width, fig_height]); % [left, bottom, width, height] 
        exportgraphics(gcf, f_name+'.eps', 'ContentType', 'vector')
        % exportgraphics(figure(idx), f_name+'.eps',"Padding","figure")
    
        % matlab2tikz(char(f_name+".tex"))
    end
end

%% QUANTITATIVE STUDY
C1_err = L2_norm(recs{3}.x_hist - xd_hist, dt);
C2_err = L2_norm(recs{2}.x_hist - xd_hist, dt);
C3_err = L2_norm(recs{1}.x_hist - xd_hist, dt);

fprintf("\nL2 Norm of Tracking Error:\n");
fprintf("C1: %.4f\n", C1_err);
fprintf("C2: %.4f\n", C2_err);
fprintf("C3: %.4f\n", C3_err);
fprintf("Improvement over C3:\n");
fprintf("C1: %.2f %%\n", (C3_err - C1_err)/C3_err*100);
fprintf("C2: %.2f %%\n", (C3_err - C2_err)/C3_err*100);

%%
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

% L2 norm calculation
function y = L2_norm(x, dt)
    y = sum(sqrt(sum(x.^2,1))*dt);
end
