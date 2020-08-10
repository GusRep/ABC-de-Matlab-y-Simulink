% Este programa muestra la imagen a color, y luego muestra cada una de sus capas en tonos de grises.
% Veremos que la escencia de la imagen puede estraerse de cualquiera de las 
% capas ya que en cada una de ellas se puede apreciar la imagen real.
% Nota: los 3 archivos de muestra deben estar en el mismo PATH
% --------------------------------------------------------------------

clear
clc
I = imread('f117.jpg');	%cargamos la imagen
subplot(2,2,1);
imshow(I)                   %muestra la imagen
title('Original')

subplot(2,2,2);
imshow(I(:,:,1))			%muestra la capa roja
title('Capa RED')

subplot(2,2,3);
imshow(I(:,:,2))			%muestra la capa verde
title('Capa GREEN')

subplot(2,2,4);
imshow(I(:,:,3))			%muestra la capa azul
title('Capa BLUE')

