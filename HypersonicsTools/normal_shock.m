function [p2, rho2, T2, u2] = normal_shock(p1, rho1, T1, u1, gas_model)
  
    gamma = gas_model.gamma(T1);
    a = sqrt(gamma*gas_model.R*T1);
    M = u1/a; 
    
    % calculate energy 
    e1 = gas_model.e(T1);
      
    % set initial guess
    xo = zeros(3,1);
    [xo(1), xo(2), xo(3)] = normal_shock_CPG_ratios(M,gamma);
%     xo(1) = 1 + 2*gamma/(1+gamma) * (M^2 - 1);
%     xo(2) = (gamma + 1)*M^2/ (2 + (gamma-1)*M^2);
%     xo(3) = xo(1)/xo(2); % = [p/p1; rho/rho1; T/T1]    
    
    % b) use fsolve
    opt = optimset('Jacobian','on','Display','off');
    x = fsolve( @(x) normal_shock_iterative(x, p1, rho1, T1, e1, u1, gas_model) , xo, opt);
    
    p2=x(1)*p1;
    rho2=x(2)*rho1;
    T2=x(3)*T1;
    u2=rho1*u1/rho2;
    
end

function [F, Jacob] = normal_shock_iterative(x, p1, rho1, T1, e1, u1, gas_model)

    % calculate the function that needs to be zero
    F=zeros(3,1);
    F(1)=   ( gas_model.e( x(3)*T1) )/e1 - 1 + (p1/(rho1*e1))*( (x(1)+1)/2 )*(1/x(2) - 1);
    F(2)=   x(1)-1+(rho1*u1^2/p1)*(1/x(2)-1);
    F(3)=   x(1)-gas_model.pressure_from_rho_T(x(2)*rho1, x(3)*T1)/p1;
    
    % calculates gas properties
    Cv=gas_model.Cv(x(3)*T1);
    kT=gas_model.kT(x(1)*p1, x(2)*rho1);
    beta=gas_model.beta(x(1)*p1, x(2)*rho1);

    % a) Jacobian
    % calculate the Jacobian
    Jacob=zeros(3,3);
    Jacob(1,1)= 0.5 * (p1/rho1/e1) * (1/x(2)-1);
    Jacob(1,2)= (p1/rho1/e1) * (1/x(2)^2) * ( (x(1)-1)/2  - x(3)*(T1*beta/p1/kT) );
    Jacob(1,3)= T1/e1*Cv; % [de/dT]_rho
    Jacob(2,1)=1;
    Jacob(2,2)=-(rho1*u1^2/p1) * (1/x(2)^2);
    Jacob(2,3)=0;
    Jacob(3,1)=1;
    Jacob(3,2)= -(1/p1) * (1/x(2)/kT) ;
    Jacob(3,3)= -(1/p1) * (beta/kT); % [dp/dT]_rho
    
end

function [p_ratio, rho_ratio, T_ratio] = normal_shock_CPG_ratios(M1,gamma)
    p_ratio = (1 + 2*gamma/(1+gamma) * (M1^2 - 1));
    rho_ratio = ((gamma + 1)*M1^2/ (2 + (gamma-1)*M1^2));
    T_ratio = p_ratio/rho_ratio;
end
