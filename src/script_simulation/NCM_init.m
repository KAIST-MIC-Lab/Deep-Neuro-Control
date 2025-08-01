
% default parameters
ncm.init = 1;       % initialization flag
ncm.dt = ctrl_dt;   % control time step
ncm.x_num = 2;      % number of state variables

ncm.alpha = 1e-5;   % decay rate (contracting)
ncm.d_MAX = 1;      % maximum disturbance (not used in this version)

% initial values
ncm.mu = 1e3;  
ncm.W_bar = 1e0*eye(ncm.x_num);  
% ncm.W_bar = ncm.mu*eye(ncm.x_num);  
ncm.X = 1e1;

% control gains
R = diag([1e0, 1e0])*1e5;
ncm.inv_R = inv(R);  

ncm.lbd = 1e-4; % penalty term for metric amplitude

% ncm.m_MAX = 0e3;    % maximum mass (not tunable)
% ncm.m_MIN = 0e-1;   % minimum mass (not tunable)
