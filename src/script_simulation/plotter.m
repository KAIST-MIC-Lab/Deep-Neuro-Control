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

%% COLORS
c_list = [
    0 0 1       % BLUE
    0 1 1       % CYAN
    .5 .5 .5    % GRAY
];

%% MAIN PLOT FUNCTIONS
figure(1); clf; 
% hF = gcf; 
% hF.Position(3:4) = [1600, 800];
tl = tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact');

% ============================================
%     STATE VARIABLE
% ============================================
for x_idx = 1:1:num_x
    nexttile;

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

    plot(sys{idx}.x_hist(1,:), sys{idx}.x_hist(2,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(sys{idx}.x_non_hist(1,:), sys{idx}.x_non_hist(2,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on
end
plot(xd_hist(1,:), xd_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on


grid on; grid minor;
xlabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$x_2$', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0
%     ylim([minVal-len*ratio maxVal+len*ratio]);
%     xlim([minVal-len*ratio maxVal+len*ratio]);
% end

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';

% ============================================
%     CONTROL INPUT
% ============================================
for u_idx = 1:1:num_u
    nexttile;

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
    ylabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
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

    plot(sys{idx}.u_hist(1,:), sys{idx}.u_hist(2,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    plot(sys{idx}.uSat_hist(1,:), sys{idx}.uSat_hist(2,:), "Color", c, "LineWidth", line_width, "LineStyle", "--"); hold on
end
plot(ud_hist(1,:), ud_hist(2,:), "Color", "red", "LineWidth", line_width, "LineStyle", "--"); hold on

rad = -pi:0.1:pi;
plot(max_u*sin(rad), max_u*cos(rad), "Color", "black")

grid on; grid minor;
xlabel('$u_1$', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('$u_2$', 'FontSize', font_size, 'Interpreter', 'latex');
% len = maxVal-minVal; ratio = .1;
% if len ~= 0
%     ylim([minVal-len*ratio maxVal+len*ratio]);
%     xlim([minVal-len*ratio maxVal+len*ratio]);
% end

ax = gca;
ax.FontSize = font_size; 
ax.FontName = 'Times New Roman';


%% 
figure(2); clf; 
% hF = gcf; 
% hF.Position(3:4) = [1600, 800];
tl = tiledlayout(4, 2, 'Padding', 'none', 'TileSpacing', 'compact');


% ============================================
%     Tracking error 
% ============================================
for x_idx = 1:1:num_x
    nexttile;

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
    ylabel('$x_1$', 'FontSize', font_size, 'Interpreter', 'latex');
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

maxVal = 0; minVal = 0;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    plot(t, sys{idx}.M_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on
    
    maxVal = max(maxVal, max(sys{idx}.M_hist(1,:)));
    minVal = min(minVal, min(sys{idx}.M_hist(1,:)));
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

maxVal = 0; minVal = 0;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    plot(t, sys{idx}.optDone_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(sys{idx}.optDone_hist(1,:)));
    minVal = min(minVal, min(sys{idx}.optDone_hist(1,:)));
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

maxVal = 0; minVal = 0;
for idx = 1:1:num_sample
    c = c_list(idx, :);

    plot(t, sys{idx}.mu_hist(1,:), "Color", c, "LineWidth", line_width, "LineStyle", "-"); hold on

    maxVal = max(maxVal, max(sys{idx}.mu_hist(1,:)));
    minVal = min(minVal, min(sys{idx}.mu_hist(1,:)));
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

%% LOCAL FUNCTIONS
function y = time_norm(x)
    num_x = size(x, 1);

    y = zeros(1,size(x,2));
    for x_idx = 1:1:num_x
        y = y + x(x_idx,:).^2;
    end
    y = sqrt(y);
end