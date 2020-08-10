function [reout,imt,w] = nyquist1(a,b,c,d,iu,w)
%NYQUIST1 Diagrama  de Nyquist de la  respuesta en frecuencia
%         para sistemas lineales de tiempo continuo.
%
%       Esta versión del comando NYQUIST tiene en cuenta los polos en
%       el eje jw y los rodea al crear el vector frecuencia  para
%       para producir el Diagrama de Nyquist correcto (El comando
%       NYQUIST no lo hace, por lo que produce un diagrama incorrecto
%       en esos casos).   Como un rasgo agregado, esta función devuelve
%       la cantidad de polos a lazo abierto en el semiplano derecho,
%       la cantidad de rodeos en contra-reloj,
%       y la cantidad de polos a lazo cerrado en el semiplano derecho
%       en su pantalla.
%
%       NOTA: Esta versión de NYQUIST1 no tiene en cuenta cancelaciones
%       polo-cero.  Por lo tanto, el ususario debe simplificar la
%       función de transferencia ante de usar este comando.
%
% USOS  :
%       NYQUIST(A,B,C,D,IU) produce un diagrama de Nyquist de
%       la entrada IU a todas las salidas del sistema:
%               .                                               -1
%               x = Ax + Bu                        G(s) = C(sI-A) B + D
%               y = Cx + Du
%                               RE(w) = real(G(jw)), IM(w) = imag(G(jw))
%
%       El rango de frecuencias y numero de puntos se elige automaticamente.
%
%       NYQUIST1(NUM,DEN) produce un diagrama de Nyquist mediante el poli-
%       nomio función de transferencia G(s) = NUM(s)/DEN(s), donde
%       NUM y DEN contienen los coeficientes en potencias descendentes de s.
%
%       NYQUIST1(A,B,C,D,IU,W) o NYQUIST(NUM,DEN,W) usa el vector
%       de frecuencias suministrado por el usuario (en rad/seg.) para
%       evaluar la respuesta en esas frecuencias.  Cuando se la invoque
%       con argumentos de salida,
%               [RE,IM,W] = NYQUIST(A,B,C,D,...)
%               [RE,IM,W] = NYQUIST(NUM,DEN,...)
%       devuelve el vector frecuencia W y las matrices RE e IM con tantas
%       columnas como salidas y tantos renglones como length(W) .
%       No grafica nada.
%       Vea también: LOGSPACE,MARGIN,BODE, y NICHOLS.

%       J.N. Little 10-11-85
%       Revised ACWG 8-15-89, CMT 7-9-90, ACWG 2-12-91, 6-21-92,
%               AFP 2-23-93
%               WCM 8-31-97
%
%  ********************************************************************
%  Modifications made to the nyquist - takes into account poles on jw axis.
%  then goes around these to make up frecuencia vector
%


if nargin==0, eval('exresp(''nyquist'')'), return, end

% --- Determina la sintaxis usada ---
if (nargin==1),
        error('número de entradas no válido.');

elseif (nargin==2),     % Función de Transferencia sin W
        num  = a; den = b;
        w = freqint2(num,den,30);
        [ny,nn] = size(num); nu = 1;

elseif (nargin==3),     % Función de Transferencia con W
        num = a; den = b;
        w = c;
        [ny,nn] = size(num); nu = 1;

% elseif (nargin==4),     % Espacio de estado , sin iu o W
%         error(abcdchk(a,b,c,d));
%         w = freqint2(a,b,c,d,30);
%         [iu,nargin,re,im]=mulresp('nyquist',a,b,c,d,w,nargout,0);
%         if ~iu, if nargout, reout = re; end, return, end
%         [ny,nu] = size(d);

elseif (nargin==5),     % Espacio de estado, con iu y sin W
        error(abcdchk(a,b,c,d));
        w = freqint2(a,b,c,d,30);
        [ny,nu] = size(d);

else
        error(abcdchk(a,b,c,d));
        [ny,nu] = size(d);

end

if nu*ny==0, im=[]; w=[]; if nargout~=0, reout=[]; end, return, end

% *********************************************************************
% disgrsión del archivo nyquist original
% tenemos un vector frecuencia, un numerador y denominador
% se crea código que rodea polos y ceros.

if (nargin==5) | (nargin ==4) | (nargin == 6)
        [num,den]=ss2tf(a,b,c,d)
end
tol = 1e-6;  % tolerancia para polos imaginarios
z = roots(num);
p = roots(den)
% ***** Si todos los polos están en el origen,
%       tan solo los movemos un poco a la izquierda***
if norm(p) == 0
 length_p = length(p);
 p = -tol*ones(length_p,1);
 den = den(1,1)*[1  tol];
 for ii = 2:length_p
  den = conv(den,[1  tol])
 end
end

zp = [z;p];        % combina los polos y ceros del sistema
nzp = length(zp);  % numero de polos y ceros
ones_zp=ones(nzp,1);
%Ipo = find((abs(real(p))<1e-6) & (imag(p)>=0));
%orden de polos con parte real nula + parte imag. no negativa
Ipo = find((abs(real(p))< tol) & (imag(p)>=0));
if  ~isempty(Ipo)   %
% **** en caso de tener dichos polos se hace lo siguiente:******
po = p(Ipo); % polos con real part 0 y parte imag. no menor que 0
% verifica los distintos tipos de polo
[y,ipo] = sort(imag(po));  % ordena por parte imaginaria
po = po(ipo);
dpo = diff(po);
idpo = find(abs(dpo)> tol);
idpo = [1;idpo+1];   % indiza los polos

po = po(idpo);   % se usan los polos distintos
nIpo = length(idpo); % # de tales polos
originflag = find(imag(po)==0);  % busca polo en el origen

s = [];  % s es el vector respuesta en frecuencia
w = sqrt(-1)*w;  % crea un vector jwo para evaluar las frecuencias

for ii=1:nIpo % todos los polos Ipo

        [nrows,ncolumns]=size(w);
        if nrows == 1
                w = w.';  % si w es una fila, lo hace columna
        end;
        if nIpo == 1
                r(ii) = .1;
        else            % distancias a otros polos y ceros
                pdiff = zp-po(ii)*ones_zp  % busca diferencias entre
                                           % el polo a verificar y los
                                           % otros polos y ceros
                ipdiff = find(abs(pdiff)> tol) % ipdiff son todas las
                                               % diferencias no nulas

                r(ii)=0.2*min(abs(pdiff(ipdiff))) % toma la mitad
                r(ii)=min(r(ii),0.1)  % y el mínimo entre la dif. y .1
        end;
        t = linspace(-pi/2,pi/2,25);
        if (ii == originflag)
                t = linspace(0,pi/2,25);
        end;    % devuelve vector de puntos alrededor de cada Ipo
        s1 = po(ii)+r(ii)*(cos(t)+sqrt(-1)*sin(t));  % desvío
        s1 = s1.';  % se asegura que es columna

        % Ahora se reconstruyen las frecuencias s complejas y se evalúa de nuevo
        bottomvalue = po(ii)-sqrt(-1)*r(ii);  % el valor de la parte imag.
        if (ii ==  originflag)  % si este es un punto de origen
        bottomvalue = 0;
        end;
        topvalue = po(ii)+sqrt(-1)*r(ii); % el último valorthe donde
                                          % concluye el desvío
        nbegin = find(imag(w) < imag(bottomvalue)); %
        nnbegin = length(nbegin); % encuentra los menores que rodean las
        if (nnbegin == 0)& (ii ~= originflag)    % raíces en el eje jw
                sbegin = 0
        else sbegin = w(nbegin);
        end;
        nend = find(imag(w) > imag(topvalue));  % puntos mayores que
        nnend = length(nend);    % rodean raíces en jw
        if (nnend == 0)
                send = 0
        else send = w(nend);
        end
        w = [sbegin; s1; send];  % reconstruye la mitad del eje jw
end
else
        w = sqrt(-1)*w;
end
%end  % esto termina el lazo para polos en el eje imaginario
% *************************************************************
% fin de disgresión, pasa al comando nyquist estándar
% Calcula la respuesta en frecuencia
if (nargin==2)|(nargin==3)
        gt = freqresp(num,den,w);
else
        gt = freqresp(a,b,c,d,iu,w);
end
% ***********************************************************

%        nw = length(gt);
%        mag = abs(gt);   % agrega factor e escala
%        ang = angle(gt);
%        mag = log2(mag+1);    % escala por log2(mag)

%        for n = 1:nw
%                h(n,1) = mag(n,1)*(cos(ang(n,1))+sqrt(-1)*sin(ang(n,1)));
%        end;  % recalcula G(jw) con el factor de escala

%        gt = h;
% ***********************************************************
ret=real(gt);
imt=imag(gt);

% Si no se usó argumentos de salida dibuja el diagrama.
if nargout==0,
   status = ishold;
   plot(ret,imt,'r-',ret,-imt,'g--')

%  plot(real(w),imag(w))
% modificaciones agregadas aquí Nueva Discgresión
        % ****************************************
          % para contar los rodeos
        [numc,denc] = tfchk(num,den);
          % crea los reflejos + y - de G(jw) en el eje imaginario
        gw = [(ret + j*imt); numc(1)/denc(1); (flipud(ret) - j*flipud(imt))]; %
          % mira G(jw)
        [Ngw,Mgw] = size(gw);  % size del G(jw) evaluado
        gwp1 = gw + ones(Ngw,Mgw);
          % mueve el origen de 0 a -1, para poder contar los rodeos a -1
        initial_angle = rem(180/pi*angle(gwp1(1)) + 360, 360);
          % define el ángulo inicial
        angle_gwp1 = rem(180/pi*angle(gwp1) - initial_angle + 720,360);
          % angulos del origen a todos los puntos
        dagw = diff(angle_gwp1);
          % resta ángulo inicial halla grados de rodeo
        tolerance = 180;
          % define tolerancia - donde el contador de rodeos "echa" afuera
        Ipd = find(dagw < -tolerance);
          % el numero de rodeos en contra-reloj
        Ind = find(dagw > tolerance);
          % este es el numero de rodeos en a-reloj
        Nacw = max(size(Ipd)) - max(size(Ind));   % el numero de rodeos
        % el criterio de Nyquist dice Z = P - N
        P = length(find(p>0))
        disp('P = números de polos a lazo abierto en el SD (rhp)');
        N = Nacw
        disp('N =  numero de rodeos en contra-reloj')
         Z = P - N     %
        disp('Z = números de polos a lazo cerrado en el SD (rhp)')
          %*******************************************


   set(gca, 'YLimMode', 'auto')
   limits = axis;
   % Pone hold on porque la próxima figura puede re-escalarse
   set(gca, 'YLimMode', 'auto')
   set(gca, 'XLimMode', 'manual')
   hold on
   % Crea flechas
   for k=1:size(gt,2),
        g = gt(:,k);
        re = ret(:,k);
        im = imt(:,k);
        sx = limits(2) - limits(1);
        [sy,muestra]=max(abs(2*im));
        arrow=[-1;0;-1] + 0.75*sqrt(-1)*[1;0;-1];
        muestra=muestra+(muestra==1);
        reim=diag(g(muestra,:));
        d=diag(g(muestra+1,:)-g(muestra-1,:));
        % Gira la flecha teneindo en cuenta los factores de escala sx y sy
        d = real(d)*sy + sqrt(-1)*imag(d)*sx;
        rot=d./abs(d);          % cuando la flecha no es horizontal
        arrow = ones(3,1)*rot'.*arrow;
        scalex = (max(real(arrow)) - min(real(arrow)))*sx/50;
        scaley = (max(imag(arrow)) - min(imag(arrow)))*sy/50;
        arrow = real(arrow)*scalex + sqrt(-1)*imag(arrow)*scaley;
        xy =ones(3,1)*reim' + arrow;
        xy2=ones(3,1)*reim' - arrow;
        [m,n]=size(g);
        hold on
        plot(real(xy),-imag(xy),'r-',real(xy2),imag(xy2),'g-')
   end
   xlabel('Real Axis'), ylabel('Imag Axis')

   limits = axis;
     % Crea cruces en  s = -1 + j0, p.e. el punto (-1,0)
   if limits(2) >= -1.5  & limits(1) <= -0.5 % Dibuja si -1 no está alejado.
        line1 = (limits(2)-limits(1))/50;
        line2 = (limits(4)-limits(3))/50;
        plot([-1+line1, -1-line1], [0,0], 'w-', [-1, -1], [line2, -line2], 'w-')
   end

   % Ejes
   plot([limits(1:2);0,0]',[0,0;limits(3:4)]','w:');

   if ~status, hold off, end    % Vuelve el hold al estado previo
   return % suprime la salida
end
%reout = ret;
%   plot(real(p),imag(p),'x',real(z),imag(z),'o');
