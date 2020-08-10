function Julia_SET_CPU(C)


%Conjunto de Julia calculado sobre el CPU. Familia Z=Z^2+C.
%El conjunto de Julia se forma haciendo las iteraciones Z=Z^2+C sobre el
%plano complejo, los puntos del plano que se van al infinito no
%representan (en Matlab esos puntos tienen la maxima amplitud, están 
%"a fondo de escala").
%C: Numero complejo. 

%fecha: 9/10/2009.
%Ejemplo: 
%Julia_SET_CPU(-1.3+0.1i).

tam=input('Rango visible del plano complejo (tipicamente 3): ');
Puntos=input('Cantidad de puntos en una dimensión del plano complejo (2000 en GPU es tipico): ');
Pasos=input('Cantidad de iteraciones (30 es lo usual): ' );


%Eje_R y Eje_Im se usan para generar dos matrices que contienen las
%coordenadas x e y del plano complejo (Acá se carga el rango del eje que se 
%quiere ver y la cantidad de puntos en el rango).
Eje_R=linspace(-tam,tam,Puntos);
Eje_Im=linspace(-tam,tam,Puntos);
%Meshgrid construye los planos de coordenadas x e y a partir de los ejes R
%e Im.
[Zr,Zi]=meshgrid(Eje_R,Eje_Im);

Zi=i.*Zi;
Z=Zr+Zi;

[m,n]=size(Zr);
Cz=C.*ones(m,n);

tic;
for j=1:1:Pasos
    Z=Z.^2+Cz;
end
disp('duración del calculo (s): ');disp(toc);


%Ajustando la amplitud para que haya sensibilidad cerca del cero.
Amp_ajustada=exp(-abs(Z));
imagesc(Amp_ajustada);
axis('image','off');

end

