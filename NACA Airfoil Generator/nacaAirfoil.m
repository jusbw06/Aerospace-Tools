function [airfoil] = nacaAirfoil(name)
%NACA MPTT
    if length(name) > 4 || length(name) < 4
       perror('Invalid NACA Airfoil') 
    end

    M = str2double(name(1))/100;
    P = str2double(name(2))/10;
    TT = str2double(name(3:4))/100;
    
    xc = 0:.01:1;
    yc = zeros(length(xc),1);
    dt = zeros(length(xc),1);
    xu = zeros(length(xc),1);
    yu = zeros(length(xc),1);
    xl = zeros(length(xc),1);
    yl = zeros(length(xc),1);
    
    for i = 1:length(xc)
        
        yc(i) = camber_eqn(M,P,xc(i));
        dt(i) = thickness_eqn(TT,xc(i));
        [xu(i), yu(i), xl(i), yl(i)] = combine_eqn(M,P,dt(i),yc(i),xc(i));
        
    end
    
    airfoil = toStruct(xu, yu, xl, yl, xc, yc dt)
    
end

% check for compatibility
function [s] = toStruct(varargin)

    for i = 1:length(varargin):
        s.{varargin{i}.name} = varargin{i}.value
    end

end


function [y] = camber_eqn(m,p,x)
        
        if x < p
            y = (m/p^2)*(2*p*x-x^2);
        elseif x >= p
            y = ( m/(1-p)^2 )*( (1-2*p) + 2*p*x - x^2);
        else
            perror('Check your variables.')
        end
        
end

function [dt] = thickness_eqn(t,x)
        
       dt = (t/0.20)*(0.2969*sqrt(x) - 0.126*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);
        
end

function [xu, yu, xl, yl] = combine_eqn(m,p,dt,yc, xc)

    m = camber_eqn_derivative(m,p,xc);
    ds = sqrt(m^2 + 1);

    temp = [xc;yc] + (dt/ds)*[-m;1];

    xu = temp(1); yu = temp(2);
    
    temp = [xc;yc] - (dt/ds)*[-m;1];

    xl = temp(1); yl = temp(2);
    
end

function [dydx] = camber_eqn_derivative(m,p,x)
        
        if x < p
            dydx = (m/p^2)*(2*p - 2*x);
        elseif x >= p
            dydx = ( m/(1-p)^2 )*( 2*p - 2*x);
        else
            perror('Check your variables.')
        end
        
end
