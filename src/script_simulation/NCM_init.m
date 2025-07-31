ncm.init = 1;

ncm.dt = ctrl_dt;
ncm.m_MAX = 0e3;  % maximum mass
ncm.m_MIN = 0e-1;   % minimum mass
ncm.d_MAX = 1;  % maximum disturbance

ncm.mu = 1e0;  
ncm.alpha = 1e-1;  % decay rate

ncm.x_num = 2;

ncm.W_bar = ncm.mu*eye(ncm.x_num);  % initial controller gain
ncm.X = ncm.mu;

R = diag([1, 1])*1e4;
ncm.inv_R = inv(R);  % input weight matrix