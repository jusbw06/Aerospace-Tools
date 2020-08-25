function [ p2, rho2, T2, M2 ] = normal_shock_CPG( p1, rho1, T1, M1, gas_model)

    gamma = gas_model.gamma(T1);
    [p_ratio, rho_ratio, T_ratio] = normal_shock_CPG_ratios(M1,gamma);
    p2 =  p_ratio * p1;
    rho2 =  rho_ratio * rho1;
    T2 = (p_ratio/rho_ratio) * T1; 
    M2 = sqrt(   (1+((gamma-1)/2)*M1^2)/(gamma*M1^2-(gamma-1)/2)   );
    
end

function [p_ratio, rho_ratio, T_ratio] = normal_shock_CPG_ratios(M1,gamma)

    p_ratio = (1 + 2*gamma/(1+gamma) * (M1^2 - 1));
    rho_ratio = ((gamma + 1)*M1^2/ (2 + (gamma-1)*M1^2));
    T_ratio = p_ratio/rho_ratio;
    
end