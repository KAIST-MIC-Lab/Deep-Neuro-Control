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
tl = tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     STATE VARIABLE
% ============================================
for x_idx = 1:1:num_x
    nexttile;

    maxVal = -9999; minVal = 9999;
    for c_idx = 1:1:case_num
        plot(t, recs{c_idx}.x_hist(x_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
        plot(t, recs{c_idx}.x_non_hist(x_idx,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "--"); hold on

        maxVal = max(maxVal, max([recs{c_idx}.x_hist(x_idx,:), recs{c_idx}.x_non_hist(x_idx,:), xd_hist(x_idx,:)]));
        minVal = min(minVal, min([recs{c_idx}.x_hist(x_idx,:), recs{c_idx}.x_non_hist(x_idx,:), xd_hist(x_idx,:)]));
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
% nexttile;
% for c_idx = 1:1:case_num
%     plot3(recs{c_idx}.x_hist(1,:), recs{c_idx}.x_hist(2,:), recs{c_idx}.x_hist(3,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
%     plot3(recs{c_idx}.x_non_hist(1,:), recs{c_idx}.x_non_hist(2,:), recs{c_idx}.x_non_hist(3,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "--"); hold on
% end
% plot3(xd_hist(1,:), xd_hist(2,:), xd_hist(3,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
% 
% grid on; grid minor;
% xlabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
% zlabel('$x_3$', 'FontSize', font_size, 'Interpreter', 'latex');
% % len = maxVal-minVal; ratio = .1;
% % if len ~= 0
% %     ylim([minVal-len*ratio maxVal+len*ratio]);
% %     xlim([minVal-len*ratio maxVal+len*ratio]);
% % end
% 
% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';
% % pbaspect([1 1 1])
% axis equal
% % axis square

% ============================================
%     CONTROL INPUT
% ============================================
for u_idx = 1:1:num_u
    nexttile;

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

% ============================================
%     CONTROL INPUT LOCI
% ============================================
% nexttile;
% max_u = param.max_u;
% 
% for c_idx = 1:1:case_num
%     plot3(recs{c_idx}.u_hist(1,:), recs{c_idx}.u_hist(2,:), recs{c_idx}.u_hist(3,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
%     plot3(recs{c_idx}.uSat_hist(1,:), recs{c_idx}.uSat_hist(2,:), recs{c_idx}.uSat_hist(3,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "--"); hold on
% end
% [x y z] = sphere(128);
% h = surfl(max_u*x, max_u*y, max_u*z); 
% set(h, 'FaceColor', [.5 .5 .5], 'FaceAlpha', 0.1, 'EdgeColor', 'none')
% % shading interp
% % plot3(ud_hist(1,:), ud_hist(2,:), ud_hist(3,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on
% 
% % rad = -pi:0.1:pi;
% % plot3(max_u*sin(rad), max_u*cos(rad), zeros(size(rad)), "Color", "black")
% 
% grid on; grid minor;
% xlabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
% ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
% zlabel('$u_3$', 'FontSize', font_size, 'Interpreter', 'latex');
% % len = maxVal-minVal; ratio = .1;
% % if len ~= 0
% %     ylim([minVal-len*ratio maxVal+len*ratio]);
% %     xlim([minVal-len*ratio maxVal+len*ratio]);
% % end
% 
% ax = gca;
% ax.FontSize = font_size; 
% ax.FontName = 'Times New Roman';
% % pbaspect([1 1 1])
% axis equal
% % axis square

%% 
figure(2); clf; 
% hF = gcf; 
% hF.Position(3:4) = [1600, 800];
tl = tiledlayout(4, 2, 'Padding', 'none', 'TileSpacing', 'compact');

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
ylabel("$\Vert u\Vert$", 'FontSize', font_size, 'Interpreter', 'latex');
maxVal = param.max_u; minVal = 0;
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';


% ============================================
%     Tracking error 
% ============================================
for x_idx = 1:1:num_x
    nexttile;

    maxVal = -9999; minVal = 9999;
    for c_idx = 1:1:case_num

        % err = time_norm(recs{c_idx}.x_hist-xd_hist);
        err = recs{c_idx}.x_hist(x_idx,:) - xd_hist(x_idx,:);
        plot(t, err, "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on

        maxVal = max(maxVal, max(err));
        minVal = min(minVal, min(err));
    end

    grid on; grid minor;
    xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
    ylabel("$e_"+x_idx+"$", 'FontSize', font_size, 'Interpreter', 'latex');
    len = maxVal-minVal; ratio = .1;
    if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
    xlim([0 T])

    ax = gca;
    ax.FontSize = font_size; 
    ax.FontName = 'Times New Roman';
end

% ============================================
%     X
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    plot(t, recs{c_idx}.M_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    
    maxVal = max(maxVal, max(recs{c_idx}.M_hist(1,:)));
    minVal = min(minVal, min(recs{c_idx}.M_hist(1,:)));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$X$', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     optDone
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    plot(t, recs{c_idx}.optDone_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    maxVal = max(maxVal, max(recs{c_idx}.optDone_hist(1,:)));
    minVal = min(minVal, min(recs{c_idx}.optDone_hist(1,:)));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Done', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     mu
% ============================================
nexttile;

maxVal = -9999; minVal = 9999;
for c_idx = 1:1:case_num
    plot(t, recs{c_idx}.mu_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    maxVal = max(maxVal, max(recs{c_idx}.mu_hist(1,:)));
    minVal = min(minVal, min(recs{c_idx}.mu_hist(1,:)));
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('mu', 'FontSize', font_size, 'Interpreter', 'latex');
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
for c_idx = 1:1:case_num
    plot(t, recs{c_idx}.contraction_flag_hist(1,:), "Color", recs{c_idx}.color, "LineWidth", line_width, "LineStyle", "-"); hold on
    maxVal = 1; minVal = 0;
end
grid on; grid minor;
xlabel('Time / s', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Contraction Condition', 'FontSize', font_size, 'Interpreter', 'latex');
len = maxVal-minVal; ratio = .1;
if len ~= 0; ylim([minVal-len*ratio maxVal+len*ratio]);  end
xlim([0 T])

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

%% LOCAL FUNCTIONS
function y = time_norm(x)
    num_x = size(x, 1);

    y = zeros(1,size(x,2));
    for x_idx = 1:1:num_x
        y = y + x(x_idx,:).^2;
    end
    y = sqrt(y);
end