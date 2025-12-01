function ncm = NCM_init(ctrl_dt, ctrl_opt)

    % default parameters
    ncm.init = 1;       % initialization flag
    ncm.dt = ctrl_dt;   % control time step
    ncm.x_num = 3;      % number of state variables

    % initial values
    mu = 1e-3;  
    W_bar = eye(ncm.x_num);  
    % W_bar = mu*eye(ncm.x_num);  
    M = inv(W_bar) * mu; % W_bar / mu = W
    % W_bar = zeros(ncm.x_num);
    % M = eye(ncm.x_num);
    X = 1e-9;

    % control gains
    R = diag([1e0 1e0 1e0])*1e0;
    ncm.inv_R = inv(R);  

    switch ctrl_opt
        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 11 
            ncm.ctrl_no = 1;
                
            ncm.alpha = 1e-3;    % decay rate (contracting)

            ncm.M = M;
            % ncm.M_bar = eye(ncmj.x_num);
            ncm.mu = mu;
            ncm.X = X;
            ncm.y = ncm.M*ones(ncm.x_num,1);

        % PROPOSED 2 – inverse space (not maintained)
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 21 % proposed 2 (inverse)
            error("not maintained anymore")

            % ncm.ctrl_no = 2;
            % ncm.AUX = eye(ncm.x_num);

        % EXISTING – penalty method
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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


    ncm.u = zeros(ncm.x_num,1);
end