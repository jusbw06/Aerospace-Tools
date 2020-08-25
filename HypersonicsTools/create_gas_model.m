function gas_model = create_gas_model(model_type, gas_type)

    if gas_type == 0 % Air
        
        gamma = 1.4;
        R = 287;
        Tv = 3521;
        
    elseif gas_type == 1 % Nitrogen
        
        gamma = 1.4;
        R = 296.8;
        Tv = 3521;
        
    else
        disp('Need more gas types!!!')
    end
    gas_model.type = model_type;


    
  if (model_type==0) % CPG
      
    Cv = R/(gamma-1);
    Cp = gamma*R/(gamma-1);
    gas_model.R = R;
    gas_model.Cv = @(T) Cv;
    gas_model.Cp = @(T) Cp;
    gas_model.gamma = @(T) gamma;
    gas_model.e = @(T) Cv*T;
    gas_model.h = @(T) Cp*T;
    gas_model.a_from_h = @(h) sqrt( gamma*R*h/Cp );
    gas_model.T_from_h = @(h) h/Cp; 
    gas_model.pressure_from_rho_T = @(rho,T)  rho*R*T;
    gas_model.kT = @(p, rho) 1/p;
    gas_model.beta = @(p, rho) rho*R/p;
    gas_model.a = @(T) sqrt(gas_model.gamma(T)*gas_model.R*T);
    
    
  elseif (model_type == 1) % TPG
      
    gas_model.R = R;
    psi = @(T) gas_model.R * (Tv/2/T / (sinh(Tv/2/T)));
    gas_model.Cv = @(T) 5/2*gas_model.R + psi(T) ;
    gas_model.Cp = @(T) gas_model.Cv(T) + gas_model.R; %7/2*gas_model.R + psi(T);
    gas_model.gamma = @(T) gas_model.Cp(T)/gas_model.Cv(T);
    gas_model.e = @(T) 5/2 * gas_model.R * T + (R*Tv)/(exp(Tv/T)-1);
    gas_model.h = @(T) gas_model.e(T) + R*T;
    gas_model.a_from_h = @(h) sqrt( gas_model.gamma( T_from_h(h,gas_model) ) * gas_model.R * T_from_h(h, gas_model) );
    gas_model.a = @(T) sqrt(gas_model.gamma(T)*gas_model.R*T);
    gas_model.T_from_h = @(h) T_from_h(h, gas_model); 
    gas_model.kT = @(p, rho) 1/p;
    gas_model.beta = @(p, rho) rho*R/p;
    gas_model.pressure_from_rho_T = @(rho,T)  rho*R*T;    
    
  else
    error('I need a new model!!!!');
  end

end


function [T] = T_from_h(h, gas_model)
    func = @(T) gas_model.h(T) - h;
    T = fzero(func, 500);
end
