clc; clear; close all;

%% Problem 3

gas_model = create_gas_model(1,0);
% use T to calc lower bound 967
M = 1:20;
for i=1:length(M)
   v(i) = v_from_M(M(i),gas_model);
end

figure;
hold all;
plot(M,v*180/pi)
title('Prandtl-Meyer angler versus Mach Number - A')
xlabel('Mach Number')
ylabel('Prandtl-Meyer Angle (deg)')


function [F1, T] = v_from_M(M_t,gas_model)
    
    T = 2500;
    
    h1 = gas_model.h(T);
    
    du = 1;
    u1 = sqrt(gas_model.gamma(T)*gas_model.R*T);
    u = u1;
    a = u1;
    M = 1;
    F1 = 0;
    
    while true
        
        if M >= M_t
            return;
        end        
        
        F1 = F1 + sqrt(  (1/a^2) - (1/u^2)  )*du ;      
        
        u = u + du;
        
        h = h1 + u1^2/2 - u^2/2;
        
        T = gas_model.T_from_h(h);

        a = sqrt(gas_model.gamma(T)*gas_model.R*T);
        
        M = u/a;
                
    end

end

function [M] = M_from_v(v,gas_model)

    func = @(M) v_from_M(M,gas_model) - v;
    M = fzero(func, [1 20]);
    
end