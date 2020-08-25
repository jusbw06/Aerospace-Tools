function [p2, rho2, T2, M2] = oblique_shock(p1, rho1, T1, M1, theta, gas_model, shock_type) 
% 1 for weak
% 0 for strong

    wave_angle = beta_from_theta(theta,M1,shock_type);
    Mn1 = sin(wave_angle)*M1;
    [p2, rho2, T2, Mn2] = normal_shock_CPG(p1, rho1, T1, Mn1, gas_model);
    M2 = Mn2/sin(wave_angle-theta);


end

function [beta] = beta_from_theta(theta, M1, delta)
% In with rad out with rad
% Do not enter 0 for theta
% delta = 1 for weak shock, delta = 0 for strong shock

    gamma = 1.4;
    lambda = ( (M1^2 - 1)^2 - 3*(1 + (gamma-1)/2*M1^2)*(1+(gamma+1)/2*M1^2)*tan(theta)^2)^(1/2); 
    x = ( (M1^2 - 1)^3 - 9*(1 + (gamma-1)/2*M1^2 )*( 1 + (gamma-1)/2*M1^2 + (gamma+1)/4*M1^4)*tan(theta)^2)/lambda^3;
    beta = atan( ( M1^2 - 1 + 2*lambda*cos( ( 4*pi*delta + acos(x) )/3 ) ) / ( 3*(1+ (gamma-1)/2 *M1^2) * tan(theta) ) ) ;

end