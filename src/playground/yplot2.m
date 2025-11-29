clear
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));
addpath(genpath([pwd filesep 'mosek/11.0/tools/platform/osxaarch64/bin']));
addpath(genpath([pwd filesep 'mosek/11.0/toolbox/r2022bom']));
% C:\Program Files\Mosek\11.0\tools\platform\win64x86\bin

%%

%% VISUALIZATION (y-space)
figure(1);clf; 
e = [1;.5];

% contraction constraint
H = -eye(2);
f = e'*[
    -1 0
    0 1
] * 1e2;
g = -1;
% H = -H; f = -f; g = -g;
fp = fimplicit(@(x,y) H(1,1)*x.^2 + (H(1,2)+H(2,1))*x.*y + H(2,2)*y.^2 + f(1)*x + f(2)*y + g, 'Color', 'black', 'LineWidth', 1); hold on
% fill(fp.XData, fp.YData, 'blue')

grid on

maxX = max(fp.XData); minX = min(fp.XData);
maxY = max(fp.YData); minY = min(fp.YData);
xlim([minX-0.1*(maxX-minX), maxX+0.1*(maxX-minX)]);
ylim([minY-0.1*(maxY-minY), maxY+0.1*(maxY-minY)]);

plot([-1e9 1e9], [-1e9 1e9]); hold on
