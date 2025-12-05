function ncm = NCM_init(ctrl_dt, ctrl_opt)

    % default parameters
    ncm.init = 1;       % initialization flag
    ncm.dt = ctrl_dt;   % control time step
    ncm.x_num = 3;      % number of state variables

    % initial values
    nu = 1e-3;  
    W_bar = eye(ncm.x_num);  
    M = (W_bar\eye(size(W_bar))) * nu; % W_bar / nu = W
    X = 1;

    % control gains
    R = diag([1e0 1e0 1e0])*1e0;
    ncm.inv_R = R\eye(size(R));

    switch ctrl_opt
        % NOMINAL
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 0
            ncm.ctrl_no = 0;
            ncm.M = 0; ncm.X = 1;

        % PROPOSED 1 – effective space
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 11 
            ncm.ctrl_no = 1;
                
            ncm.alpha = 1e1;    % decay rate (contracting)

            ncm.M_energy = 3e1;
            ncm.M = eye(ncm.x_num)*ncm.M_energy;

            ncm.nu = 0;
            ncm.X = X;
            ncm.y = ncm.M*ones(ncm.x_num,1);

            ncm.lbd = 50e1 ; % penalty term
            % ncm.lbd = 50e-2 ; % penalty term
        
            % EXISTING – CV-STEM
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        case 31 % large penalty -> poor tracking
            ncm.ctrl_no = 3;

            ncm.alpha = 1e1;    % decay rate (contracting)

            ncm.W_bar = W_bar;
            ncm.nu = nu;
            ncm.M = M;    
            ncm.X = X;

            ncm.lbd = 1e-1 ; % penalty term
        
        case 32 % small penalty -> better but satuatrion occurs
            ncm.ctrl_no = 3;

            ncm.alpha = 1e1;    % decay rate (contracting)

            ncm.W_bar = W_bar;
            ncm.nu = nu;   
            ncm.M = M;    
            ncm.X = X;

            ncm.lbd = 1e-10; % penalty term
        
        otherwise
            error("Invalid control option")
    end


    ncm.u = zeros(ncm.x_num,1);
end