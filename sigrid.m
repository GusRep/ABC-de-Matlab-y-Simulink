function[ ] = sigrid(sig)

%SIGRID  Genera l�neas de trama en el plano s para rlocus o pzmap.
%
%       SIGRID genera una grilla en el plano s sobre  un Lugar de ra�ces
%       o mapa polo-zero existente. Se dibujan las l�neas a sigma constante.
%       �sese con SGRID si sigma, zeta, y Wn se requieren en forma
%       simult�nea.  Puede usarse sola si se desea.
%
%       Vea Tambi�n: RLOCUS, ZGRID, SGRID, JGRID, y PZMAP.

error(nargchk(1,1,nargin));

hold on

%Plot sigma line
    limits = axis;
    mx=limits(1,4);
    mn=limits(1,3);
    stz=abs(mx)+abs(mn);
    st=stz/50;
    im=mn:st:mx;
    lim=length(im);
    for i=1:lim
        re(i)=-sig;
    end
    re(:);
    plot(re,im,'.')

hold off

return
