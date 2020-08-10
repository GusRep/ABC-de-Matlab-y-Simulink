
% ejemplo de filtrado “Smooth” empleando la foto de Lenna

G = imread('lenna.jpg'); 		%cargamos la imagen

% Matriz máscara que determinará el efecto buscado
HSM=[1/16 8/16 1/16; 8/16 16/16 8/16; 1/16 8/16 1/16]; 

A=rgb2gray(G); 		% Transforma en escala de grises
A=double(A);   		% Transforma a A en "double"

C=conv2(A,HSM); 		% Filtro Smooth
A=uint8(A); 			% Transforma a A en “uint8”
C=uint8(C); 			% Transforma a C en “uint8”

subplot (1,2,1);

imshow (A); 			% Imagen base sin aplicar filtro
title ('Imagen sin filtrar');


subplot (1,2,2);
imshow (C); 			% Imagen con filtro "Smooth"
title ('Imagen con filtro "Smooth"');
