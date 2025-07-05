ncm.dt = dt;
ncm.m_MAX = 1e1;  % maximum mass
ncm.m_MIN = 1e0;   % minimum mass
ncm.d_MAX = 1;  % maximum disturbance

ncm.mu = 1e-1;  
ncm.alpha = 1e1;  % decay rate
ncm.W = zeros(4);  % initial controller gain

R = 1e1;
ncm.inv_R = inv(R);  % input weight matrix