% Demostración "interactiva" en Matlab de EFICIENCIA COMPUTACIONAL
echo off
% Inicializamos
clear
clc
format 
% Nada de lo que escribamos con ECHO OFF puede verse, LALALA !!!
echo on
% Tema: EFICIENCIA COMPUTACIONAL by Gus
% Demostración "interactiva" en Matlab
%
% ¿COMO AUMENTAR LA VELOCIDAD Y MEJORAR EL MANEJO DE MEMORIA?
%  =========================================================
% 
% Las operaciones vectoriales y matriciales en MATLAB, son mas
% rápidas que las que se hacen elemento a elemento.  
%
% Esto quiere decir que para obtener la máxima velocidad es 
% conveniente vectorizar los algoritmos en los archivos ".m".
%
% Por tanto es deseable, siempre que sea posible, expresar los
% bucles "for" y "while", mediante operaciones vectoriales o matriciales.
%
% 
% EJEMPLOS:
% =========
% A continuación escribimos dos secuencias distintas que permiten
% calcular el seno de 1001 números.
% Comprendidos entre 0 y 10.
% 
%	i=0;
%	for t=0:.01:10          <--- Definimos un vector de 0 a 10 con
%	                             incrementos de 0.01 (los 1001)
%		i=i+1;
%		y(i)=sin(t);
%	end
% 
% y una versión vectorizada de lo mismo es:
% 
%	t=0:.01:10;
%	y=sin(t);               <--- Definimos un vector "y" de resultados
% 
% Pulsar una tecla para continuar...
pause
clc

% 	Vamos a calcular los tiempos de cálculo de CPU de una y otra
% rutina de código, para lo cual conviene hacer dos archivos  
% con el siguiente contenido:
% 
clear
format long
t0=cputime;
i=0;
for t=0:.01:10
	i=i+1;          
	y(i)=sin(t);
end            
t1=cputime-t0;
disp(['Tiempo de cpu (versión secuencial): ' num2str(t1)])
% 
%   Pulsar una tecla para continuar...
pause
% =====================================================================
%	Y para la versión vectorizada:
% 
clear
format long
t0=cputime;
t=0:.01:10;
y=sin(t)           % <--- Luego correlo si el  ; Que debría haber?
t1=cputime-t0;
disp(['Tiempo de cpu (version vectorizada): ' num2str(t1)])
% 
%	Con la versión secuencial el tiempo de cálculo, en un 486DX2-66MHz,
% en el mejor de los casos es de 0.44 sg., mientras que con la versión
% vectorizada, en el peor de los casos es de 0.05.  Como puede apreciarse
% la diferencia en tiempo de ejecución es ABISMAL !!! 
% (del orden de 10 veces menor con la versión vectorizada).
% 
% Analizar lo que sucedió en nuestra PC.
%
% Teniendo presenete que no muestra en pantalla el proceso inermedio.
%
% Pulsar una tecla para continuar...
pause
clc

% VECTORES PREDEFINIDOS
% =====================
%	Si no se puede vectorizar una parte del código, se puede conseguir
% que los bucles "for" se ejecuten mas rápidamente, simplemente
% predefiniendo vectores en los cuales se almacenen los datos de salida.
% Por ejemplo, se puede incluir como primera sentencia una función "zeros",
% con lo cual el bucle "for" se ejecuta un poco mas rápido

clear
format long
x=[1 2 3 4;5 6 3 4;6 -3 4 2;5 7 9 1];
max=5000;        % Cambiando el valor de max hacemos mas o menos cálculos
t0=cputime;
y=zeros(1,max);
for i=1:max
	y(i)=det(x^i);
end
t1=cputime-t0;
disp(['Tiempo de cpu. ' num2str(t1)])

%   Pulsar una tecla para continuar...
pause
clc

%   que con la misma versión sin predefinir el vector y

clear
format long
x=[1 2 3 4;5 6 3 4;6 -3 4 2;5 7 9 1];
max=5000;        % Cambiando el valor de max hacemos mas o menos cálculos
t0=cputime;
for i=1:max
	y(i)=det(x^i);
end
t1=cputime-t0;
disp(['Tiempo de cpu. ' num2str(t1)])

% CONCLUSION: si no se predefine el vector y, MATLAB debe redimensionar 
% dicho vector a un elemento mas en cada iteración (algo que no concierne 
% a nuestro programa, es algo no vital).
% En otras palabras, decirle previamente al MATLAB, la cantidad  
% de memoria que voy a necesitar.
%
% Si se predefine el vector, se elimina esta etapa "no vital" y se 
% ejecuta mas rápido.
%
% FIN
echo off