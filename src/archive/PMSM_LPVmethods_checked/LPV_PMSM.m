%% LPV method for PMSM
%% LPV methods [2023 - 05 - 18]

Ae0 = [0 1 0 0;0 -B/J 0 0;0 0 -R/L 0;0 0 0 -R/L];
Be = [0 0 1/L 0;0 0 0 1/L]';
Ahate1 = [0 0 0 0;0 0 -1/J 0;0 1/L 0 0;0 0 0 0];
Ahate2 = [0 0 0 0; 0 0 0 1/J;0 0 0 0;0 -1/L 0 0];

%% Polytopic decomposed matrices: Expanded

eta1_low = -Km;     % nonlienar term lower boundary
eta1_high = Km;     % nonlienar term higher boundary
eta2_low = -Km;     % nonlienar term lower boundary
eta2_high = Km;     % nonlienar term higher boundary


% Ae1 = Ae0 + eta1_low * Ahate1;
% Ae2 = Ae0 + eta1_high * Ahate1;
% Ae3 = Ae0 + eta2_low * Ahate2;
% Ae4 = Ae0 + eta2_high * Ahate2;


% Expansion coefficient: 2
Ae1 = Ae0 + 2 * eta1_low * Ahate1;          
Ae2 = Ae0 + 2 * eta1_high * Ahate1;
Ae3 = Ae0 + 2 * eta2_low * Ahate2;
Ae4 = Ae0 + 2 * eta2_high * Ahate2;

sigma1 = 0.5;
sigma2 = 0.5;


% Vhat = [sqrt(2)*[eta1_low eta1_high 0 0;0 0 eta2_low eta2_high];1 1 0 0;0 0 1 1];
% Vhat = [[eta1_low eta1_high 0 0;0 0 eta2_low eta2_high];0.5 0.5 0 0;0 0 0.5 0.5];
Vhat = [[eta1_low eta1_high 0 0;0 0 eta2_low eta2_high];1 1 0 0;0 0 1 1];
Vhatinv = inv(Vhat);

syspole = [-40 -80 -100 -120];
K1 = place(Ae1,Be,syspole);
K2 = place(Ae2,Be,syspole);
K3 = place(Ae3,Be,syspole);
K4 = place(Ae4,Be,syspole);