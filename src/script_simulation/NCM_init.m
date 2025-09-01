function ncm = NCM_init(ctrl_dt, ctrl_opt)

    % default parameters
    ncm.init = 1;       % initialization flag
    ncm.dt = ctrl_dt;   % control time step
    ncm.x_num = 3;      % number of state variables

    ncm.alpha = 1e0;   % decay rate (contracting)
    ncm.d_MAX = 1;      % maximum disturbance (not used in this version)

    % initial values
    ncm.mu = 1e0;  
    ncm.W_bar = 1e0*eye(ncm.x_num);  
    % ncm.W_bar = ncm.mu*eye(ncm.x_num);  
    ncm.X = 1e1;

    % control gains
    R = diag([1e0 1e0 1e0])*1e0;
    ncm.inv_R = inv(R);  

    switch ctrl_opt
        case 1 % proposed
            ncm.cstr_on = 1;

            ncm.AUX = eye(ncm.x_num);
            
        case 2 % large penalty -> poor tracking
            ncm.cstr_on = 0;

            ncm.lbd = 1e-0 ; % penalty term
        
        case 3 % small penalty -> better but satuatrion occurs
            ncm.cstr_on = 0;

            ncm.lbd = 1e-3; % penalty term
        
        otherwise
            error("Invalid control option")
    end


    % ncm.m_MAX = 0e3;    % maximum mass (not tunable)
    % ncm.m_MIN = 0e-1;   % minimum mass (not tunable)

end