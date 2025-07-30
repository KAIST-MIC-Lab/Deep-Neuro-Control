
function SDC = SDCmotor(x,xd)

    %%
    % L = .66e-3;    % Inductance (mH)
    L = .66;
    R = 0.251;     % Resistance (Ohm)
    J = 3.24e-5;   % Inertia (kg.m^2)
    Phi = 16.8e-3; % Flux (Wb)
    P = 4;         % Pole pairs
    fv = 2e-3;     % Viscous friction (N.m.s/rad)

    %%
    a21 = 1/J * (-(3/2)*P*Phi) * ...
        (-intXcosY(x(3),xd(3), x(1),xd(1), P) + ...
        intXsinY(x(4),xd(4), x(1),xd(1), P));
    a23 = 1/J * (-(3/2)*P*Phi) * ...
        intXsinY(1,1, x(1),xd(1), P);
    a24 = 1/J * (+(3/2)*P*Phi) * ...
        intXcosY(1,1, x(1),xd(1), P);

    a31 = 1/L*P^2*Phi * ...
        intXcosY(x(2),xd(2), x(1),xd(1), P);
    a41 = 1/L*P^2*Phi * ...
        intXsinY(x(2),xd(2), x(1),xd(1), P);

    a32 = 1/L*P*Phi * ...
        intXsinY(1,1, x(1),xd(1), P);
    a42 = 1/L*P*Phi * ...
        intXcosY(1,1, x(1),xd(1), P);

    %%
    SDC = [
        0 1 0 0;
        a21 -fv/J a23 a24;
        a31 a32 -R/L 0;
        a41 a42 0 -R/L
    ];

    for i = 1:1:4
        for j = 1:1:4
            if isnan(SDC(i,j))
                SDC(i,j) = 0;
            end
        end
    end
end

function z = intXsinY(x1,x2, y1,y2, a)
    z = 1/(a*(y1-y2)) * (-cos(a*y1)*x1 + cos(a*y2)*x2) + (x1-x2)/(a*(y1-y2))^2 * (sin(a*y1) - sin(a*y2));
end

function z = intXcosY(x1,x2, y1,y2, a)
    z = 1/(a*(y1-y2)) * (sin(a*y1)*x1 - sin(a*y2)*x2) + (x1-x2)/(a*(y1-y2))^2 * (-cos(a*y1) + cos(a*y2));
end