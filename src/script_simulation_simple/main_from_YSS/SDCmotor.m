
function SDC = SDCmotor(x,xd,ud,param)
    x1 = x(1);  % theta
    % x2 = x(2);  % theta dot

    r1 = xd(1);  % reference position
    % dr1 = xd(2);  % reference velocity
    
    ud1 = ud(1);  % control input 1
    ud2 = ud(2);  % control input 2

    %%
    % L = param.L;    % Inductance (mH)
    P = param.P;         % Pole pairs
    Phi = param.Phi; % Flux (Wb)
    J = param.J;   % Inertia (kg.m^2)
    fv = param.fv;     % Viscous friction (N.m.s/rad)

    %%
    if x1 ~= r1
        b2 = (3*P*Phi)/(2*J*(x1-r1)) * ...
            ( ...
                - (sin(P*x1) - sin(P*r1)) * ud1 ...
                + (cos(P*x1) - cos(P*r1)) * ud2 ...
            );
    else
        b2 = (3*P*Phi)/(2*J) * ...
            ( ...
                - P*(cos(P*x1)) * ud1 ...
                - P*(sin(P*x1)) * ud2 ...
            );
    end

    %%
    SDC = [
        0 1;
        b2 -fv/J
    ];

    % for i = 1:1:2
    %     for j = 1:1:2
    %         if isnan(SDC(i,j))
    %             SDC(i,j) = 0;
    %         end
    %     end
    % end
end
