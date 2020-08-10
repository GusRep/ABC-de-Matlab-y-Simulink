function [reout,imt,w] = lnyquist1(a,b,c,d,iu,w)
%LNYQUIST1 Nyquist respuesta en frecuencia para sistemas lineales de tiempo continuo.
%
%       Esta Versión del  Comando LNYQUIST tiene en cuenta polos en el
%       eje jw y los rodea cuando crea  el vector frecuencia  para
%       para producir el diagrama de Nyquist (el comando NYQUIST no
%       lo hace y por lo tanto produce un diagrama incorrecto cuando tenemos 
%       polos en el eje jw).
%       Como característica extra, esta función devuelve la cantidad de
%       polos a lazo abierto en el semiplano derecho ,
%       la cantidad de rodeos en contra-reloj, y la cantidad de
%       polos a lazo cerrado en el semiplano derecho en su pantalla.
%
%       NOTA: Esta versión de LNYQUIST1 no tiene cuenta las cancelaciones 
%       polo-cero.  Por lo tanto, el usuario debe simplificar la función
%       de transferencia antese de usar este comando.
%      

%LNYQUIST1 respuesta en frecuencia para sistemas lineales de tiempo continuo.
%       NYQUIST(A,B,C,D,IU) produce un diagrama de Nyquist de 
%       la entrada IU a todas las salidas del sistema:
%               .                                    -1
%               x = Ax + Bu             G(s) = C(sI-A) B + D
%               y = Cx + Du      RE(w) = real(G(jw)), IM(w) = imag(G(jw))
%
%       El rango de frecuencias y número de puntos se eligen automaticamente.
%
%       LNYQUIST1(NUM,DEN) produce el diagrama de Nyquist para la función de
%       transferencia polinomial G(s) = NUM(s)/DEN(s) donde NUM y DEN con-
%       tiene los coeficientes del polinomio en potencias descendentes de s.
%
%       LNYQUIST1(A,B,C,D,IU,W) o LNYQUIST(NUM,DEN,W) usa el vector freq. W
%       suministrado por el usuario el cual debe contener las frecuencias, en 
%       radianes/seg, a las cuales se va a evaluar la respuesta.  
%       Cuando es invocado con argumentos a izquierda,
%               [RE,IM,W] = LNYQUIST(A,B,C,D,...)
%               [RE,IM,W] = LNYQUIST(NUM,DEN,...)
%       devuelve el vector frecuencia W y las matrices RE y IM con tantas
%       columnas como salidas y length(W) renglones.  No se dibuja nada en 
%       pantalla.
%       Vea además: LOGSPACE,MARGIN,BODE, y NICHOLS.

%       J.N. Little 10-11-85
%       Revised ACWG 8-15-89, CMT 7-9-90, ACWG 2-12-91, 6-21-92,
%               AFP 2-23-93
%               WCM 8-30-97
%
%  ********************************************************************
%  Modificaciones a nyquist - tiene en cuenta polos en eje jw.
%  entonces los rodea para formar el vector frecuencia
%


if nargin==0, eval('exresp(''nyquist'')'), return, end

% --- Determina la sintaxis usada ---
if (nargin==1),
        error('número de argumentos de entrada no válido.');

elseif (nargin==2),     % Función de Transferencia sin vector frecuencia
        num  = a; den = b;
        w = freqint2(num,den,30);
        [ny,nn] = size(num); nu = 1;

elseif (nargin==3),     % Función de Transferencia con vector frecuencia
        num = a; den = b;
        w = c;
        [ny,nn] = size(num); nu = 1;

% elseif (nargin==4),     % Espacio de estado , sin iu o vector frecuencia
%         error(abcdchk(a,b,c,d));
%         w = freqint2(a,b,c,d,30);
%         [iu,nargin,re,im]=mulresp('nyquist',a,b,c,d,w,nargout,0);
%         if ~iu, if nargout, reout = re; end, return, end
%         [ny,nu] = size(d);

elseif (nargin==5),     % Espacio de estado , con iu pero sin v.freq. 
        error(abcdchk(a,b,c,d));
        w = freqint2(a,b,c,d,30);
        [ny,nu] = size(d);

else
        error(abcdchk(a,b,c,d));
        [ny,nu] = size(d);

end

if nu*ny==0, im=[]; w=[]; if nargout~=0, reout=[]; end, return, end

% *********************************************************************
% Aquí se produce una disgrsión respecto del original
% tenemos un vector frecuencia,  numerador y denominador
% creamos código para rodear polos y ceros acá.

if (nargin==5) | (nargin ==4) | (nargin == 6)
        [num,den]=ss2tf(a,b,c,d)
end
tol = 1e-6;  % tolerancia para encontrar polos imaginarios
z = roots(num);
p = roots(den);
% ***** Si todos los polos están en el origen, 
%       tan solo los movemos un cachito a la izquierda***
if norm(p) == 0
 length_p = length(p);
 p = -tol*ones(length_p,1);
 den = den(1,1)*[1 tol];
 for ii = 2:length_p
  den = conv(den,[1 tol])
 end
end
zp = [z;p];        % combine los polos y ceros del sistema
nzp = length(zp);  % número de polos y ceros
ones_zp=ones(nzp,1);
%Ipo = find((abs(real(p))<1e-6) & (imag(p)>=0)) %indiza polos con cero real part + non-neg parte imag.
Ipo = find((abs(real(p)) < tol) & (imag(p)>=0)); %indiza polos con cero real part + non-neg parte imag.
if  ~isempty(Ipo)   %
% **** caso si tenemos tales polos hacemos lo siguiente:*************************
po = p(Ipo); % polos con parte real 0 y parte imag. no negativa
% check for distinct polos
[y,ipo] = sort(imag(po));  % ordena por parte imaginaria
po = po(ipo);
dpo = diff(po);
idpo = find(abs(dpo)>tol);
idpo = [1;idpo+1];   % indexes del distinct polos

po = po(idpo);   % solo se usan los polos distintos
nIpo = length(idpo); % # de tales polos
originflag = find(imag(po)==0);  % busca polo en origen

s = [];  % s es la respuesta en vector frecuencia
w = sqrt(-1)*w;  % create a jwo vector to evaluate all frecuencias with
for ii=1:nIpo % for all Ipo polos

  [nrows,ncolumns]=size(w);
  if nrows == 1
          w = w.';  % si w es un renglón, lo hace columna
  end;
  if nIpo == 1
          r(ii) = .1;
  else            % ver. distancias a los otros polos y ceros
          pdiff = zp-po(ii)*ones_zp  % halla diferencias entre
          % polos en verif. y otros polos y ceros
          ipdiff = find(abs(pdiff)>tol) % ipdiff todas dif. no nulas

          r(ii)=0.2*min(abs(pdiff(ipdiff))) % la mitad de eso
          r(ii)=min(r(ii),0.1)  % el mínimo entre esta dif. y 0.1
  end;
  t = linspace(-pi/2,pi/2,25);
  if (ii == originflag)
          t = linspace(0,pi/2,25);
  end;    % nos da un vector de puntos alrededor de cada Ipo
  s1 = po(ii)+r(ii)*(cos(t)+sqrt(-1)*sin(t));  % rodeo
  s1 = s1.';  % asegura es columna

  % reconstituye s - frecuencia compleja - 
  % y recalcula.

  bottomvalue = po(ii)-sqrt(-1)*r(ii);  % magnitud de parte imag.
  if (ii ==  originflag)  % si es un origen 
          bottomvalue = 0;
  end;
  topvalue = po(ii)+sqrt(-1)*r(ii); % ultimo valor rodeo se detiene
  nbegin = find(imag(w) < imag(bottomvalue)); %
  nnbegin = length(nbegin); % puntos menores de rodeos
  if (nnbegin == 0)& (ii ~= originflag)    % a polos en jw 
          sbegin = 0
  else sbegin = w(nbegin);
  end;
  nend = find(imag(w) > imag(topvalue));  % puntos mayores de rodeos
  nnend = length(nend);    % a raices en jw
  if (nnend == 0)
          send = 0
  else send = w(nend);
  end
  w = [sbegin; s1; send];  % medio eje jw reconstruido
end
else
  w = sqrt(-1)*w;
end
%end  % this ends el lazo for eje imaginario polos
% *************************************************************
% fin de disgresión
% calcula respuesta en frecuencia
if (nargin==2)|(nargin==3)
        gt = freqresp(num,den,w);
else
        gt = freqresp(a,b,c,d,iu,w);
end
% ***********************************************************

        nw = length(gt);
        mag = abs(gt);   % factor de escala
        ang = angle(gt);
        mag = log2(mag+1);    % escalado en log2(mag) 

        for n = 1:nw
                h(n,1) = mag(n,1)*(cos(ang(n,1))+sqrt(-1)*sin(ang(n,1)));
        end;  % recalcula G(jw) con factor de escala

        gt = h;
% ***********************************************************
ret=real(gt);
imt=imag(gt);

% Si ha llamado sin argumentos a izquierda , dibuja.
if nargout==0,
   status = ishold;
   plot(ret,imt,'r-',ret,-imt,'g--')

%  plot(real(w),imag(w))
% aquí : modificaciones adicionadas

   % ****************************************
     % aquí : modificaciones adicionadas para contar rodeos
   [numc,denc] = tfchk(num,den);
     % crea  + y - reflección de G(jw) en eje imag.
   gw = [(ret + j*imt); numc(1)/denc(1); (flipud(ret) - j*flipud(imt))]; %
     %mira G(jw)
   [Ngw,Mgw] = size(gw);  % size de G(jw) filas y cols
   gwp1 = gw + ones(Ngw,Mgw);
   % mueve el origen desde 0 a -1, para poder
     %  contar rodeos a -1 
     initial_angle = rem(180/pi*angle(gwp1(1)) + 360, 360);
     % define ángulo inicial
   angle_gwp1 = rem(180/pi*angle(gwp1) - initial_angle + 720,360);
     % ángulo del origen a los puntos
   dagw = diff(angle_gwp1);
     % resta ángulo inicial p/hallar grados de rodeo
   tolerance = 180;
     % define tolerancia - donde contador de rodeos "se da vuelta" 
   Ipd = find(dagw < -tolerance);
     % número de rodeos en contra-reloj
   Ind = find(dagw > tolerance);
     % número de rodeos a-reloj
   Nacw = max(size(Ipd)) - max(size(Ind));   % número de rodeos
   % Criterio Nyquist  Z = P - N
   P = length(find(p>0))
   disp('P = número de polos a lazo abierto en rhp (spd)');
   N = Nacw
   disp('N =  número de rodeos en contra-reloj')
    Z = P - N     %
   disp('Z = número de polos a lazo cerrado en rhp')
     %*******************************************


 set(gca, 'YLimMode', 'auto')
 limits = axis;
 % mantiene ejes porque proximo gráfico podrá reescalarlo
 set(gca, 'YLimMode', 'auto')
 set(gca, 'XLimMode', 'manual')
 hold on
 % Flechitas
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
    % gira flecha t. en cuenta factores sx y sy
    d = real(d)*sy + sqrt(-1)*imag(d)*sx;
    rot=d./abs(d);          % cuando flecha no horizontal
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
 xlabel('Eje Real'), ylabel('Eje Imag.')

 limits = axis;
 % cruz en s = -1 + j0,(punto -1)
 if limits(2) >= -1.5  & limits(1) <= -0.5 % dibuja solo si -1 point no lejos.
      line1 = (limits(2)-limits(1))/50;
      line2 = (limits(4)-limits(3))/50;
      plot([-1+line1, -1-line1], [0,0], 'w-', [-1, -1], [line2, -line2], 'w-')
 end

 % Ejes
 plot([limits(1:2);0,0]',[0,0;limits(3:4)]','w:');

 if ~status, hold off, end    % vuelve hold al estado anterior
 return % no muestra la salida
end
% reout = ret;
%   plot(real(p),imag(p),'x',real(z),imag(z),'o');
