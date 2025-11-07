% MATLAB Code for Plotting the Exterior of an Ellipsoid
% YALMIP 툴박스와 SDPT3/MOSEK 같은 솔버가 설치되어 있어야 합니다.

yalmip('clear');
clear;
clc;

%% 1. 문제 정의 (행렬 A, H와 벡터 x 정의)
% 예시로 2x2 행렬과 2x1 벡터를 사용합니다.
% u가 2차원 평면에 그려집니다.
n = 2;
A = [0.5 1; -0.5 0.2];
H = [2 0.5; 0.5 1.5]; % H는 반드시 양의 정부호 행렬이어야 함
x = [1; 1];

%% 2. 최적화 변수 설정
M = sdpvar(n, n, 'full');

%% 3. 경계선 추적을 위한 최적화
thetas = linspace(0, 2*pi, 100); % 0도부터 360도까지 탐색
boundary_points = []; % 경계점들을 저장할 배열

fprintf('경계점을 계산하는 중...\n');
for i = 1:length(thetas)
    w = [cos(thetas(i)); sin(thetas(i))]; % 탐색 방향 벡터
    
    % 제약 조건: '금지된 영역'의 경계 또는 내부
    % -M'HM + AM + M'A' >= 0
    Constraints = [-M'*H*M + A*M + M'*A' >= 0];
    
    % 목적 함수: w 방향으로 가장 멀리 있는 점 찾기 (최대화)
    Objective = -w' * (M*x); % YALMIP은 최소화를 기본으로 하므로 음수를 붙여 최대화
    
    % 최적화 실행
    ops = sdpsettings('solver', 'mosek', 'verbose', 0); % MOSEK 또는 'sdpt3' 사용
    sol = optimize(Constraints, Objective, ops);
    
    % 결과 저장
    if sol.problem == 0 % 해를 찾았다면
        M_sol = value(M);
        u_sol = M_sol * x;
        boundary_points = [boundary_points, u_sol];
    end
end
fprintf('계산 완료.\n');

%% 4. 결과 시각화
figure;
hold on;
grid on;
axis equal;

% 계산된 경계점들로 '금지된 영역' (타원체)을 채워서 그리기
fill(boundary_points(1,:), boundary_points(2,:), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1.5);

% 그래프 범위 설정 및 꾸미기
xlims = [min(boundary_points(1,:))-2, max(boundary_points(1,:))+2];
ylims = [min(boundary_points(2,:))-2, max(boundary_points(2,:))+2];
xlim(xlims);
ylim(ylims);

title('u의 해 공간 (회색 영역)');
xlabel('u_1');
ylabel('u_2');
legend('금지된 영역');

% 회색으로 바깥 영역(해 공간)을 표시
patch([xlims(1) xlims(2) xlims(2) xlims(1)], [ylims(1) ylims(1) ylims(2) ylims(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.5);
% 다시 금지된 영역을 위로 올림
fill(boundary_points(1,:), boundary_points(2,:), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'LineWidth', 1.5);

hold off;