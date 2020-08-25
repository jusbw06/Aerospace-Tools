function [p2, rho2, T2, M2] = expansion_shock(p1, rho1, T1, M1, theta, gas_model)
% CPG only


    gamma = gas_model.gamma(T1);
    
    v1 = v_from_M( M1, gamma );
    v2 = v1 + theta;
    M2 = M_from_v( v2, gamma );

    p_r1 = ( 1+(gamma-1)/2*M1^2 )^(gamma/(gamma-1));
    rho_r1 = ( 1+(gamma-1)/2*M1^2 )^(1/(gamma-1));
    T_r1 = 1 + (gamma-1)/2*M1^2;
    
    p_r2 = ( 1+(gamma-1)/2*M2^2 )^(gamma/(gamma-1));
    rho_r2 = ( 1+(gamma-1)/2*M2^2 )^(1/(gamma-1));
    T_r2 = 1 + (gamma-1)/2*M2^2;
    
    T2 = T1*T_r1*(1/T_r2);
    p2 = p1*p_r1*(1/p_r2);
    rho2 = rho1*rho_r1*(1/rho_r2);

    
end

function [PM_angle] = v_from_M(M,gamma)
% outputs rad
% only for CPG

    PM_angle = sqrt( (gamma + 1) / (gamma - 1) ) * atan( sqrt( (gamma - 1) / (gamma + 1) * (M^2 - 1) ) ) - atan( sqrt( M^2 - 1 ) );

end

function [M] = M_from_v(v,gamma)

    func = @(M) sqrt( (gamma + 1) / (gamma - 1) ) * atan( sqrt( (gamma - 1) / (gamma + 1) * (M^2 - 1) ) ) - atan( sqrt( M^2 - 1 ) ) - v;
    M = fzero(func, 3);
    
end
