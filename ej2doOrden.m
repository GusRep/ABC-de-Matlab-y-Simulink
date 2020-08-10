% Ejemplo: ej2doOrden.m
%
% En este ejemplo variaremos el parametro zeta en un sist.de 2ºorden
% ----------------------------------------------------------------------------

clear all;

% Frecuencia de oscilación natural
wn=1;
% Variaremos zeta (Coef. de amortiguamiento)
zeta=[0 0.3 0.5 1 2];
% Tiempo de análisis
t=0:20/201:20;

% Por cada parametro zeta almacenaremos la respuesta
% del sistema en una matriz y:
for n=1:length(zeta)
    num=wn^2;
    den=[1 2*zeta(n)*wn wn^2];
    W=tf(num,den);
    y(:,n)=step(W,t);
end;

%Representacion en 2D de la respuesta
figure(1)
plot(t,y,'LineWidth',1.8);xlabel('t');ylabel('y(t)');
title('Respuesta de un sistema de segundo orden variando \zeta');
legend('\zeta=0 (sin amortiguamiento)',...
'\zeta=0.3 (Subamortiguado)','\zeta=0.5 (Subamortiguado)',...
'\zeta=1 (Criticamente Amortiguado)','\zeta=2 (Sobreamortiguado)');
