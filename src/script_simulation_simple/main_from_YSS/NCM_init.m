% 
% % default parameters
% ncm.init = 1;       % initialization flag
% ncm.dt = ctrl_dt;   % control time step
% ncm.x_num = 2;      % number of state variables
% 
% ncm.alpha = 1e-5;   % decay rate (contracting)
% ncm.d_MAX = 1;      % maximum disturbance (not used in this version)
% 
% % initial values
% ncm.mu = 1e3;  
% ncm.W_bar = 1e0*eye(ncm.x_num);  
% % ncm.W_bar = ncm.mu*eye(ncm.x_num);  
% ncm.X = 1e1;
% 
% % control gains
% R = diag([1e0, 1e0])*1e5;
% ncm.inv_R = inv(R);  
% 
% ncm.lbd = 1e-4; % penalty term for metric amplitude
% 
% % ncm.m_MAX = 0e3;    % maximum mass (not tunable)
% % ncm.m_MIN = 0e-1;   % minimum mass (not tunable)

%% NCM controller: stabilized defaults (recommended)
ncm.init   = 1;               % initialization flag
ncm.dt     = ctrl_dt;         % control step
ncm.x_num  = 2;               % number of states

% --- Contraction / scaling ---
ncm.alpha  = 5e-3;            % contraction rate (조금 더 강하게)
ncm.mu     = 1e1;             % scale for W_bar/mu (너무 크거나 작지 않게)
ncm.mu_min = 1e-4;            % guard for runtime use
ncm.mu_max = 1e2;             % (옵션) 상한: 너무 커지는 것 방지 (필요 없으면 Inf)

% --- Metric initialization (well-conditioned PD) ---
ncm.W_bar  = 1.0*eye(ncm.x_num);  % PD metric seed
ncm.X      = 1.0;                 % gain seed (필요 시 사용)
ncm.X_min  = 1.0;                 % W_bar의 하한(I 스케일)
ncm.X_max  = 1e3;                 % (옵션) 상한: 메트릭 폭주 방지 (필요 없으면 Inf)

% --- Numerical guards ---
ncm.eps_reg  = 1e-8;   % 제어법 계산 시 W_bar+eps*I 정규화
ncm.eps_psd0 = 1e-8;   % 1차 시도 LMI PSD 여유
ncm.eps_psd1 = 1e-6;   % 실패 시 재시도용 더 큰 여유
ncm.rho_slack = 1e4;   % LMI 슬랙 s 패널티 (크게 잡아 s 사용 최소화)

% --- Input weighting (R too large -> 수치 나쁨) ---
R = 1e3*eye(2);                 % 이전 1e5 -> 1e3로 완화
ncm.inv_R = diag(1./diag(R));   % inv() 대신 안전한 방식

% --- Metric magnitude penalty ---
ncm.lbd = 1e-3;                 % 이전 1e-4 -> 1e-3로 살짝 강화

% (현재 버전 미사용) 물리 파라미터 가드 예시
% ncm.m_MAX = 0;
% ncm.m_MIN = 0;

% 더 강한 수축, 더 큰 mu 페널티, 메트릭 상한 타이트
ncm.alpha   = 5e-3;     % 1e-3 -> 5e-3 (수축 강화)
ncm.lbd     = 1e-2;     % 1e-3 -> 1e-2 (mu 페널티 강화)

ncm.mu_min  = 1e-3;     % 1e-4 -> 1e-3 (너무 작은 mu 회피)
ncm.mu_max  = 10;       % Inf -> 10   (mu 상한 강제)

ncm.X_min   = 0.5;      % 1.0 -> 0.5  (하한 약간 완화)
ncm.X_max   = 100;      % 1e3 -> 100  (메트릭 폭주 방지)

% 수치 여유 확대 (LMI 안정)
ncm.eps_psd0 = 1e-7;    % 1e-8 -> 1e-7
ncm.eps_psd1 = 1e-5;    % 1e-6 -> 1e-5
ncm.rho_slack= 1e4;     % 유지

% 입력 가중은 너무 크게 두면 inv_R가 작아져 mu를 더 키우게 됨 → 조금 낮춤
R = 5e3*eye(2);                 % 기존 1e3 → 5e3 (보상항이 약해짐)
ncm.inv_R = diag(1./diag(R));

% --- NCM correction gain (0~1) ---
ncm.kc = 0.3;   

%%%%%250813
ncm.alpha   = 5e-3;     % 수축 강화
ncm.lbd     = 2e-2;     % mu 페널티 강화

ncm.mu_min  = 1e-3;
ncm.mu_max  = 5;        % ★ μ 상한
ncm.X_min   = 0.5;
ncm.X_max   = 30;       % ★ 메트릭 상한

ncm.eps_psd0 = 1e-7;
ncm.eps_psd1 = 1e-5;
ncm.rho_slack= 1e4;

% 입력 가중을 크게 → 보상항 작게
R = 0.1e4*eye(2);         % ★ 2e4
ncm.inv_R = diag(1./diag(R));

% 보상 항 스케일
ncm.kc = 0.15;          % ★ 0.15 (0.1~0.3 범위 추천)
%% 

%%best20250813
ncm.mu_min  = 1e-2;     % 1e-3 → 1e-2  (너무 작은 μ 회피: K 급변 방지)
ncm.mu_max  = 4;        % 5 → 4
ncm.X_max   = 30;       % 50 → 30       (메트릭 상한 더 타이트)
ncm.alpha   = 1e-2;     % 1e-2 → 8e-3   (너무 강하지 않게 한 클릭 내림)

% 입력 가중/보상
R = 0.1e4*eye(2);         % 2e4 → 3e4 (보상 더 약하게)
ncm.inv_R = diag(1./diag(R));
ncm.kc    = 0.12;       % 0.15 → 0.12 (기본 보상 세기 낮춤)
%% 
%%best20250814
ncm.mu_min  = 1e-2;     % 1e-3 → 1e-2  (너무 작은 μ 회피: K 급변 방지)
ncm.mu_max  = 4;        % 5 → 4
ncm.X_max   = 30;       % 30       (메트릭 상한 더 타이트)
ncm.alpha   = 0.8e-2;     % 1e-2 → 8e-3   (너무 강하지 않게 한 클릭 내림)

% 입력 가중/보상
R = 0.06e4*eye(2);         % 2e4 → 3e4 (보상 더 약하게)
ncm.inv_R = diag(1./diag(R));
