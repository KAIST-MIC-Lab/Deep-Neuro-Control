function ncm = NCM_init(ctrl_dt, ctrl_opt)

    % default parameters
    ncm.init = 1;       % initialization flag
    ncm.dt = ctrl_dt;   % control time step
    ncm.x_num = 3;      % number of state variables

    % ncm.alpha = 1e5;    % decay rate (contracting)
    ncm.d_MAX = 1;      % maximum disturbance (not used in this version)

    % initial values
    mu = 1e0;  
    % W_bar = 1e0*eye(ncm.x_num);  
    W_bar = mu*eye(ncm.x_num);  
    M = inv(W_bar) * mu;
    X = 1e1;

    % control gains
    R = diag([1e0 1e0 1e0])*1e0;
    ncm.inv_R = inv(R);  

    switch ctrl_opt
        case 11 % proposed 1 (effective space)
            ncm.ctrl_no = 1;
                
            ncm.alpha = 1e10;    % decay rate (contracting)


            ncm.M = M;
            ncm.M_bar = eye(ncm.x_num);
            ncm.mu = mu;
            ncm.X = X;

        case 21 % proposed 2 (inverse)
            error("not maintained anymore")

            % ncm.ctrl_no = 2;
            % ncm.AUX = eye(ncm.x_num);

        case 31 % large penalty -> poor tracking
            ncm.ctrl_no = 0;

            ncm.W_bar = W_bar;
            ncm.mu = mu;
            ncm.M = M;    
            ncm.X = X;

            ncm.lbd = 1e-0 ; % penalty term
        
        case 32 % small penalty -> better but satuatrion occurs
            ncm.ctrl_no = 0;

            ncm.alpha = 1e0;    % decay rate (contracting)

            ncm.W_bar = W_bar;
            ncm.mu = mu;   
            ncm.M = M;        
            ncm.X = X;

            ncm.lbd = 1e-3; % penalty term
        
        otherwise
            error("Invalid control option")
    end

end