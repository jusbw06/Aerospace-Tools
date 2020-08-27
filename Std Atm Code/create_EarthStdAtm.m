function [EarthStdAtm] = create_EarthStdAtm()
  EarthStdAtm.gamma = 1.4; %Ratio of Specific Heats
  EarthStdAtm.a0 = 340.294; %m/s Speed of Sound
  EarthStdAtm.g0 = 9.80665; %m/s^2 Acc due to Gravity
  EarthStdAtm.m0 = 28.96442; %kg/kmol Molar Mass
  EarthStdAtm.p0 = 101325; %N/m^2 Atmospheric pressure
  EarthStdAtm.R_c = 8314.32; %Nm/kmolK Universal Gas Constant
  EarthStdAtm.R = 287.05287; %Nm/kgK Gas Constant for Air
  EarthStdAtm.r_earth = 6.356766E6; %m Earth Radius
  EarthStdAtm.S = 110.4; %K Sutherland coefficient
  EarthStdAtm.T0 = 288.15; %K Temp
  EarthStdAtm.beta_s = 1.458E-6; %Ns/m^2K^.5
  EarthStdAtm.mu0 = 17.894E-6; %Ns/m^2 Dynamic Viscosity
  EarthStdAtm.rho0 = 1.225; %kg/m^2 Density
  EarthStdAtm.layersHp = [0 11000 20000 32000 47000 50000];
  EarthStdAtm.A = [288.15 216.65 196.65 139.05 270.65 ];
  EarthStdAtm.B = [-6.5E-3 0 10^-3 2.8E-3 0 ];
  EarthStdAtm.C = [8.9619638 NaN .70551848 .34926867 NaN ];
  EarthStdAtm.D = [-.20216125E-3 NaN 3.5876861E-6 7.033098E-6 NaN ];
  EarthStdAtm.E = [5.2558797 NaN -34.163218 -12.201149 NaN ];
  EarthStdAtm.F = [NaN 128244.5 NaN NaN 41828.420 ];
  EarthStdAtm.G = [NaN -.15768852E-3 NaN NaN -.12622656E-3];
  EarthStdAtm.I = [1.048840 NaN .9726309 .84392929 NaN ];
  EarthStdAtm.J = [-23.659414E-6 NaN 4.946E-6 16.993902E-6 NaN ];
  EarthStdAtm.L = [4.2258797 NaN -35.163218 -13.201149 NaN ];
  EarthStdAtm.M = [NaN 2.0621400 NaN NaN .53839563 ];
  EarthStdAtm.N = [NaN -.15768852E-3 NaN NaN -.12622656E-3];
  EarthStdAtm.osphere = @(Hp,ind) EarthStdAtm_osphere(Hp,ind,EarthStdAtm);
  EarthStdAtm.opause = @(Hp,ind) EarthStdAtm_opause(Hp,ind,EarthStdAtm);
  EarthStdAtm.indcheck = @(ind) EarthStdAtm_indcheck(ind,EarthStdAtm);
  EarthStdAtm.calc = @(Hp_vec) StdAtm_calc(Hp_vec,EarthStdAtm);
  EarthStdAtm.cas = @(tas,Hp) EarthStdAtm_cas(tas,Hp,EarthStdAtm);
  EarthStdAtm.rho_offstd = @(Hp,delT)EarthStdAtm_rho_offstd(Hp,delT,EarthStdAtm);
  
  EarthStdAtm.cas_m_Hp_plot = @(M_vec,cas_vec,Hp_vec) cas_m_Hp_plot(M_vec,cas_vec,Hp_vec,EarthStdAtm);
  %Verification of Model
  % Hp_vec = linspace(0,EarthStdAtm.layersHp(length(EarthStdAtm.layersHp)),100);
  % solution = EarthStdAtm.calc(Hp_vec);
  %
  %
  % subplot(1,3,1)
  % plot(solution.T,Hp_vec)
  % xlabel('Temperature (K)')
  % subplot(1,3,2)
  % plot(solution.P,Hp_vec)
  % title('Model Verification')
  % xlabel('Pressure (Pa)')
  % subplot(1,3,3)
  % plot(solution.rho_std,Hp_vec)
  % xlabel('Density (kg/m^3)')
  function [layer_handle] = EarthStdAtm_indcheck(ind,EarthStdAtm)
    if isnan(EarthStdAtm.F(ind)) %==1
      layer_handle = EarthStdAtm.osphere;
    else
      layer_handle = EarthStdAtm.opause;
    end
  end
  function [T,P,rho_std,a] = EarthStdAtm_osphere(Hp,ind,EarthStdAtm)
    A = EarthStdAtm.A(ind);
    B = EarthStdAtm.B(ind);
    C = EarthStdAtm.C(ind);
    D = EarthStdAtm.D(ind);
    E = EarthStdAtm.E(ind);
    I = EarthStdAtm.I(ind);
    J = EarthStdAtm.J(ind);
    L = EarthStdAtm.L(ind);
    T = A+B*Hp;
    P = (C+D*Hp)^E;
    rho_std = (I+J*Hp)^L;
    a = sqrt(EarthStdAtm.gamma*EarthStdAtm.R*T);
  end
  function [T,P,rho_std,a] = EarthStdAtm_opause(Hp,ind,EarthStdAtm)
    A = EarthStdAtm.A(ind);
    B = EarthStdAtm.B(ind);
    F = EarthStdAtm.F(ind);
    G = EarthStdAtm.G(ind);
    M = EarthStdAtm.M(ind);
    N = EarthStdAtm.N(ind);
    T = A+B*Hp;
    P = F*exp(G*Hp);
    rho_std = M*exp(N*Hp);
    a = sqrt(EarthStdAtm.gamma*EarthStdAtm.R*T);
  end
  function [rho_offstd] = EarthStdAtm_rho_offstd(Hp,delT,EarthStdAtm)
    [std] = EarthStdAtm.calc(Hp);
    rho_offstd = std.rho_std/(1+(delT/std.T));
  end
  
  function [cas] = EarthStdAtm_cas(tas,Hp,EarthStdAtm)
    static = EarthStdAtm.calc(Hp);
    M = tas/static.a;
    % if tas<static.a
    % q = .5*static.rho_std*tas^2;
    % else 
    q = .5*static.rho_std*(tas^2)*(1+(M^2)/4);
    %end
    cas = sqrt((2*EarthStdAtm.gamma/...
    (EarthStdAtm.gamma-1))*(EarthStdAtm.p0/EarthStdAtm.rho0)*(((q/...
    EarthStdAtm.p0)+1)^((EarthStdAtm.gamma-1)/EarthStdAtm.gamma) -1));
  end
  
  
  end