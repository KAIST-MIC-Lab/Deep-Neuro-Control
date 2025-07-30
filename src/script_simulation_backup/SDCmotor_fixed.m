
function SDC = SDCmotor_fixed(x,xd)

    th = x(1);
    w = x(2);
    ia = x(3);
    ib = x(4);

    %%
    % L = .66e-3;    % Inductance (mH)
    L = .66;
    R = 0.251;     % Resistance (Ohm)
    J = 3.24e-5;   % Inertia (kg.m^2)
    Phi = 16.8e-3; % Flux (Wb)
    P = 4;         % Pole pairs
    fv = 2e-3;     % Viscous friction (N.m.s/rad)

    %%
    a21 = 1/J*(-3/2*P*Phi*cos(P*th)*P*ia -3/2*P*Phi*sin(P*th)*P*ib);
    a23 = 1/J*(-3/2*P*Phi*sin(P*th));
    a24 = 1/J*(+3/2*P*Phi*cos(P*th));

    a31 = 1/L*(+P*Phi*w*cos(P*th)*P);
    a32 = 1/L*(+P*Phi*sin(P*th));

    a41 = 1/L*(-P*Phi*w*sin(P*th)*P);
    a42 = 1/L*(-P*Phi*cos(P*th));

    %%
    SDC = [
        0       1       0       0;
        a21     fv/J    a23     a24;
        a31     a32     -R/L    0;
        a41     a42     0       -R/L;  
    ];

end

