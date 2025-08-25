% =========================================================================
% LMI와 Sector Condition을 이용한 비선형 시스템 안정성 분석
% =========================================================================
% 필요 툴박스: YALMIP, SDP Solver (SeDuMi 또는 SDPT3)
% =========================================================================

clear; clc; close all;
addpath(genpath([pwd filesep 'YALMIP-master']));
addpath(genpath([pwd filesep 'sedumi']));
addpath(genpath([pwd filesep 'mosek/11.0/tools/platform/osxaarch64/bin']));
addpath(genpath([pwd filesep 'mosek/11.0/toolbox/r2022bom']));
%% 1. 시스템 정의 (Linear Part)
% x_dot = Ax + Bw
% z = Cx + Dw
A = [0 1; -1 -0.5];
B = [0; 1];
C = [1 0]; % z = x1
D = 0;

% 시스템의 개루프(Open-loop) 특성 확인 (참고용)
fprintf('>> 시스템의 개루프 고유값 (불안정 확인):\n');
disp(eig(A));

%% 2. 비선형성 정의 및 섹터 조건
% 비선형 함수 phi(z)가 섹터 [0, K]에 속한다고 가정
K = 0.8; % Sector [0, 0.8]

% 시뮬레이션에 사용할 비선형 함수 정의 (예: Saturation)
phi = @(z) min(K, max(0, z)); % 0과 K 사이에서 saturation

fprintf('\n>> 비선형성은 Sector [0, %.2f] 조건을 만족한다고 가정합니다.\n', K);

%% 3. LMI 변수 및 제약조건 정의 (YALMIP)

% LMI 변수 정의
n = size(A, 1); % 시스템 차원
P = sdpvar(n, n, 'symmetric'); % P는 대칭 행렬
tau = sdpvar(1, 1);             % tau는 스칼라

% LMI 제약조건 설정
LMI1 = (P >= 1e-9 * eye(n)); % P > 0
LMI2 = ([A'*P + P*A, P*B - tau*K*C';
         B'*P - tau*K*C, -2*tau] <= -1e-9 * eye(n+1)); % 메인 LMI < 0
LMI3 = (tau >= 0); % S-procedure 변수 조건

Constraints = [LMI1, LMI2, LMI3];

%% 4. LMI 문제 풀이
options = sdpsettings('solver', 'sedumi', 'verbose', 1);
sol = optimize(Constraints, [], options); % 목표함수 없이 feasibility 문제로 풂

%% 5. 결과 분석
if sol.problem == 0
    fprintf('\n>> LMI가 풀렸습니다 (Feasible). 시스템은 안정합니다. 🥳\n');
    % 찾은 P와 tau 값 확인
    P_sol = value(P)
    tau_sol = value(tau)
else
    fprintf('\n>> LMI를 풀 수 없습니다 (Infeasible). 안정성을 보장할 수 없습니다. 😥\n');
    disp(sol.info);
    return; % 안정하지 않으면 시뮬레이션 불필요
end

%% 6. 시뮬레이션으로 안정성 검증
fprintf('\n>> 시뮬레이션을 통해 상태 변수가 0으로 수렴하는지 확인합니다...\n');

% 비선형 폐루프 시스템 ODE 정의
% x_dot = Ax + Bw = Ax - B*phi(Cx)
ode_func = @(t, x) A*x - B*phi(C*x);

% 초기 조건 및 시뮬레이션 시간
x0 = [2; -1];
t_span = [0 30];

% ODE 풀이
[t, x] = ode45(ode_func, t_span, x0);

% 결과 시각화
figure(1); clf;
plot(t, x(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(t, x(:,2), 'r--', 'LineWidth', 1.5);
grid on;
title('상태 변수 시뮬레이션 결과');
xlabel('시간 (t)');
ylabel('상태 값 (x)');
legend('x_1', 'x_2');
set(gcf, 'color', 'w');

figure(2); clf;
plot(t, phi(C*x'), 'k-', 'LineWidth', 1.5);
grid on;
title('비선형 함수 시뮬레이션 결과');
xlabel('시간 (t)');
ylabel('비선형 함수 값 (phi)');
set(gcf, 'color', 'w');

% 제어입력
figure(3); clf;
plot(t, C*x', 'm-', 'LineWidth', 1.5);
grid on;
title('제어 입력 시뮬레이션 결과');
xlabel('시간 (t)');
ylabel('제어 입력 값 (u)');
set(gcf, 'color', 'w');
