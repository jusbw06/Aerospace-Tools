clc; clear; close all;
%% Problem 5

gas_model = create_gas_model(0,0);

alpha = 1:26;
for i=1:length(alpha)
   phi(i) = super_sonic_plate(deg2rad(alpha(i)), gas_model); 
end

figure;
plot(alpha,rad2deg(phi))
xlabel('Angle of Attack (deg)')
ylabel('\phi (deg)')
title('\phi vs \alpha')


function [phi] = super_sonic_plate(alpha, gas_model)

    M1 = 3.0;
    T1 = 300; % K
    rho1 = 1.225; % kg/m^3
    p1 = 101325;

    % Stage 1
    [pu1, rhou1, Tu1, Mu1] = expansion_shock(p1, rho1, T1, M1, alpha, gas_model); % upper
    [pl1, rhol1, Tl1, Ml1] = oblique_shock(p1, rho1, T1, M1, alpha, gas_model, 1); % lower

    % Stage 2
    phi = 0;
    pl2 = 12;
    pu2 = -20;
    while abs(abs(pu2/pl2) - 1) > 0.001
        [pu2] = oblique_shock(pu1, rhou1, Tu1, Mu1, alpha + phi, gas_model, 1); % upper
        [pl2] = expansion_shock(pl1, rhol1, Tl1, Ml1, alpha + phi, gas_model); % lower
         phi = phi + 0.0001;
    end

end