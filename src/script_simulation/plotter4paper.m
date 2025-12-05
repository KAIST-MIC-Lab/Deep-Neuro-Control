

%% FIGURE SETTING
FIGURE_SAVE_FLAG = 1;   % save the figure as .png and .eps

font_size = 20;
line_width = 2;
lgd_size = 2;

x_list = ["x", "y", "z"];
x_vec = "\xi";

%% MAIN PLOT

%% MAIN PLOT FUNCTIONS
figure(1); clf; 
hF = gcf; 
hF.Position(3:4) = [800, 500];
tl = tiledlayout(3, 2, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     STATE VARIABLE
% ============================================
for x_idx = 1:1:num_x
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
    nexttile(2*(u_idx));

    maxVal = -9999; minVal = 9999;
    for c_idx = 1:1:case_num
        plot(t, recs{c_idx}.u_hist(u_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
        plot(t, recs{c_idx}.uSat_hist(u_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "--"); hold on

        maxVal = max(maxVal, max([recs{c_idx}.u_hist(u_idx,:), recs{c_idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
        minVal = min(minVal, min([recs{c_idx}.u_hist(u_idx,:), recs{c_idx}.uSat_hist(u_idx,:), ud_hist(u_idx,:)]));
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

%% 
figure(2); clf; 
hF = gcf; 
hF.Position(3:4) = [800, 150];
tl = tiledlayout(1, 2, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     Control norm
% ============================================
nexttile;

plot([-10, T+10], param.max_u * [1 1], "Color", "black", "LineWidth", line_width, "LineStyle", "-"); hold on
maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    u_norm = time_norm(recs{c_idx}.u_hist);
    plot(t, u_norm, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(u_norm));
    minVal = min(minVal, min(u_norm));
end

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
% nexttile;
% 
% maxVal = -9999; minVal = 9999;
% for c_idx = 1:1:case_num
% 
%     err = time_norm(recs{c_idx}.x_hist-xd_hist);
%     % err = recs{c_idx}.x_hist(x_idx,:) - xd_hist(x_idx,:);
%     plot(t, err, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
% 
%     maxVal = max(maxVal, max(err));
%     minVal = min(minVal, min(err));
% end
% 
% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel("$\Vert {\mbox{\boldmath $e$}}\Vert$", 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
% xlim([0 T])
% 
% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';

% ============================================
%     X
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    plot(t, recs{c_idx}.X_hist, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    
    maxVal = max(maxVal, max(recs{c_idx}.X_hist));
    minVal = min(minVal, min(recs{c_idx}.X_hist));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$\chi$ / ', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     contraction condition satisfaction
% ============================================
% nexttile;
% 
% maxVal = -9999; minVal = 9999;
% for c_idx = 1:1:case_num
%     plot(t, recs{c_idx}.contraction_flag_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
%     maxVal = 1; minVal = 0;
% end
% grid on; grid minor;
% xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('Contraction Cond.', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
% xlim([0 T])
% 
% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';


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