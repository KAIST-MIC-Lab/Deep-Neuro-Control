ncm.init = 1;

ncm.dt = ctrl_dt;
ncm.m_MAX = 0e3;  % maximum mass
ncm.m_MIN = 0e-1;   % minimum mass
ncm.d_MAX = 1;  % maximum disturbance

ncm.mu = 1e0;  
ncm.alpha = 1e-1;  % decay rate

ncm.x_num = 3;

ncm.W = zeros(ncm.x_num);  % initial controller gain

R = diag([1, 1])*1e-3;
ncm.inv_R = inv(R);  % input weight matrix