%% 1. 시스템 및 제어기 파라미터 정의
clear; clc; close all;

% --- 일반적인 '느린' 시스템 정의 (2차 불안정 시스템) ---
A = [0 1; 2 1]; % 불안정한 시스템
B = [0; 1];
C = [1 0];

% --- 제어 파라미터 ---
u_max = 5;      % 제어 입력 포화 한계
ref = 1;        % 목표 추종 값

% --- PI 제어기 설계 (Augmented System Pole Placement) ---
% 시스템을 [x; x_i] 형태로 확장하여 PI 제어기 이득 [K, Ki] 계산
Aa = [A, zeros(2,1); -C, 0];
Ba = [B; 0];
p = [-2, -3, -4]; % 확장 시스템의 목표 극점 (안정하고 적절한 속도)
Ka = place(Aa, Ba, p);
K = Ka(1:2);
Ki = Ka(3);

% --- 안티와인드업(Anti-windup) 파라미터 ---
% 적분기를 얼마나 빠르게 되돌릴지 결정하는 이득
Kaw = 10; % Back-calculation gain

%% 2. 시뮬레이션 설정
t_span = [0 10];
x0 = [0.5; 0];          % 시스템 초기 상태
xi0 = 0;                % 적분기 초기 상태

%% 3. 시뮬레이션 실행

% --- 제어기 1: 안티와인드업(AW)이 적용된 PI 제어기 ---
y0_aw = [x0; xi0];
[t_aw, y_aw] = ode45(@(t,y) general_pi_dynamics(t,y,A,B,C,K,Ki,u_max,Kaw,ref,'AW'), t_span, y0_aw);

% --- 제어기 2: 일반 PI 제어기 (Naive) ---
y0_naive = [x0; xi0];
[t_naive, y_naive] = ode45(@(t,y) general_pi_dynamics(t,y,A,B,C,K,Ki,u_max,Kaw,ref,'Naive'), t_span, y0_naive);


%% 4. 결과 분석 및 그래프 출력
[u_aw, ua_aw] = calculate_pi_signals(t_aw, y_aw, K, Ki, u_max);
[u_naive, ua_naive] = calculate_pi_signals(t_naive, y_naive, K, Ki, u_max);

figure('Position', [100 100 1200 800]);

% 서브플롯 1: 시스템 출력 y 응답
subplot(2, 2, 1);
plot(t_aw, y_aw(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_naive, y_naive(:,1), 'r--', 'LineWidth', 2);
yline(ref, 'k:', 'LineWidth', 1.5);
grid on;
title('시스템 출력 (y = x1)');
xlabel('시간 (s)'); ylabel('y');
legend('안티와인드업 적용', '일반 제어기', '목표값', 'Location', 'SouthEast');

% 서브플롯 2: 적분기 상태 (xi)
subplot(2, 2, 2);
plot(t_aw, y_aw(:,3), 'b-', 'LineWidth', 2); hold on;
plot(t_naive, y_naive(:,3), 'r--', 'LineWidth', 2);
grid on;
title('제어기 내부 적분기 상태 (x_i)');
xlabel('시간 (s)'); ylabel('x_i');
legend('안티와인드업 적용', '일반 제어기', 'Location', 'SouthEast');

% 서브플롯 3: 제어기가 계산한 명령 (u)
subplot(2, 2, 3);
plot(t_aw, u_aw, 'b-', 'LineWidth', 2); hold on;
plot(t_naive, u_naive, 'r--', 'LineWidth', 2);
yline(u_max, 'k:', 'LineWidth', 1.5, 'Label', '포화 한계');
yline(-u_max, 'k:', 'LineWidth', 1.5);
grid on;
title('제어기가 계산한 명령 (u)');
xlabel('시간 (s)'); ylabel('u');
legend('안티와인드업 적용', '일반 제어기', 'Location', 'SouthEast');

%% 시스템 동역학 함수
function dydt = general_pi_dynamics(t, y, A, B, C, K, Ki, u_max, Kaw, ref, type)
    % 상태 변수 분리
    x = y(1:2);  % 시스템 상태
    xi = y(3);   % 적분기 상태
    
    % 제어 입력 u 계산
    u = -K * x + Ki * xi;
    
    % 입력 포화 적용
    u_a = max(min(u, u_max), -u_max);
    
    % 상태 미분 방정식
    dxdt = A * x + B * u_a;
    
    % 적분기 미분 방정식 (Anti-windup 적용 지점)
    error = ref - C * x;
    if strcmp(type, 'AW')
        % AW: 포화 오차(u-u_a)를 이용해 적분기 누적을 방지 (Back-calculation)
        dxidt = error - Kaw * (u - u_a); 
    else % Naive
        % Naive: 포화와 상관없이 오차를 계속 적분
        dxidt = error;
    end
    
    dydt = [dxdt; dxidt];
end

%% 그래프 출력을 위한 제어 신호 계산 함수
function [u, u_a] = calculate_pi_signals(t, y, K, Ki, u_max)
    u = zeros(length(t), 1);
    for i = 1:length(t)
        x = y(i, 1:2)';
        xi = y(i, 3);
        u(i) = -K * x + Ki * xi;
    end
    u_a = max(min(u, u_max), -u_max);
end