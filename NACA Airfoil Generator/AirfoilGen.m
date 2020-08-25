clc; clear; close all;

airfoil_name = '4429';

[af] = nacaAirfoil(airfoil_name);


figure;
subplot(3,1,2)
plot(af.xc,af.yc)
xlabel('Chord Length')
ylabel('Height')
title('Airfoil Camber')
axis  equal
subplot(3,1,3)
plot(af.xc,af.dt)
xlabel('Chord Length')
ylabel('Height')
title('Airfoil Thickness')
axis equal
subplot(3,1,1)
hold all;
plot(af.xu,af.yu)
plot(af.xl,af.yl, 'b')
xlabel('Chord Length')
ylabel('Height')
title(['Airfoil Dimensions -- NACA ' airfoil_name])
axis equal
    
[pass, I] = max(af.yu);
disp(['Maximum y-coordinate occurs at: x = ' num2str(af.xu(I))])
