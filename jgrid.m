function[ ] = jgrid(zeta,sigma)
%JGRID  Genera las líneas para un lugar de raíces o un mapa polo-cero
%       en el plano s
%       JGRID dibuja líneas auxiliares de razón de amortiguamiento
%	    constante (zeta) y de sigma constante en el plano continuo
%	    de un ploteo existente, ya sea obtenido por rlocus o pzmap.
%
%       Vea también: RLOCUS, ZGRID, SGRID, y PZMAP.

error(nargchk(2,2,nargin));

status = ishold;

hold on

%Plot líneas de sigma
    limits = axis;
    for i = 1:30*max(limits(:)),
        sl(i) = 1;
    end
    [w,z] = meshgrid(sigma,sl);
    re = -w.*z;
    [mcols, nrows] = size(z);
    im = 0.1:0.1:3*max(limits(:));
    plot(re,im,'w:',re,-im,'w:');
    hold on

%Plot  líneas de amortiguamiento
    [w,z] = meshgrid([0;sigma(:);2*max(limits(:))]',zeta);
    re =  -w.*z;
	[mcols, nrows] = size(z);
    im = w.*sqrt( ones(mcols,nrows) -z.*z);
	plot(re',im','w:',re',-im','w:');

if (~status), hold off, end