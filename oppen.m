% EXPERIMENTO DE OPPENHEIM
% ========================
% 1) LECTURA DE LA IMAGEN A FILTRAR
G1=imread('F117_ByN.gif');
G2=imread('Esfera_ByN.gif');
GD1=double(G1); 
GD2=double(G2);

% 2) TRANSF. DE FOURIER DE LA IMAGEN A FILTRAR
GD1=fft2(GD1);
GD2=fft2(GD2);
faseGD1=angle(GD1);
GD1abs=abs(GD1);
faseGD2=angle(GD2);
GD2abs=abs(GD2);

% 3) MEZCLA DE LOS MÓDULOS DE LAS FASES 
P1=GD1abs.*exp(i.*faseGD2);
P2=GD2abs.*exp(i.*faseGD1);

% 4) ANTITRANSF. DE FOURIER DE LA IMAGEN FILTRADA
afilt1=abs(ifft2(P1));
maxafilt1=max(max(afilt1));
H1=afilt1/maxafilt1;
afilt2=abs(ifft2(P2));
maxafilt2=max(max(afilt2));
H2=afilt2/maxafilt2;
                               
% 5) GRAFICAS
figure (1);
axis equal;
subplot (2,4,1);
imshow(G1);
title('Imagen #1: F117');
subplot (2,4,2);
imshow(GD1);
title('FFT de la imagen #1: F117');
subplot (2,4,3);
imshow(faseGD1);
title({'Fase de la FFT de';'la imagen #1: F117'});
subplot (2,4,5);
imshow(G2);
title('Imagen #2: Esfera');
subplot (2,4,6);
imshow(GD2);
title('FFT de la imagen #2: Esfera');
subplot (2,4,7);
imshow(faseGD2);
title({'Fase de la FFT de';'la imagen #2: Esfera'});
subplot (2,4,4);
imshow(H1);
title('IFFT de la imagen #1: F117');
subplot (2,4,8);
imshow(H2);
title('IFFT de la imagen #2: Esfera');
