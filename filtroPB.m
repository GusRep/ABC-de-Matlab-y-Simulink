% FILTRO PASA BAJOS
% =================
G=imread('lenna.jpg');
G=rgb2gray(G); % Este comando se ejecuta sólo si la imagen es ".jpg"
GD=double(G);
%
% 2) TRANSF. FOURIER DE LA IMAGEN
GD=fft2(GD); % Transformo la imagen
tfa=fftshift(GD);% Normalizo
tfaabs=abs(tfa); % Extraigo valores absolutos 
% maxtfaabs=max(max(tfaabs));
% PP=tfaabs/(40*sqrt(maxtfaabs));
%
% 3) VENTANA
[d1,d2]=size(GD);
ventana=zeros(size(GD));
a=98; %Base de la ventana
b=43; %Altura de la ventana
ventana((d1/2)-a:(d1/2)+a,(d2/2)-b:(d2/2)+b)=1;
%
% 4) ACCIÓN DEL FILTRO
tfafilt=tfa.*ventana;
tfafiltabs=abs(tfafilt);
% maxtfafiltabs=max(max(tfafiltabs));
% QQ=tfafiltabs/(40*sqrt(maxtfafiltabs));
%
% 5) ANTITRANSF. FOURIER DE LA IMAGEN FILTRADA
afilt=abs(ifft2(tfafilt));
maxafilt=max(max(afilt));
HH=afilt/maxafilt;
%
% 6) GRAFICA
figure;
subplot(2,3,[1 4]);
imshow(G); title('Imagen original');
%
subplot(2,3,[2 5]);
imshow(HH); title('Imagen filtrada');
%
% subplot(2,3,3);
% imhist(G);
% title('Histograma de la imagen original');
% xlabel('Intensidad de los pixeles');
% ylabel('Cantidad de pixeles');
% subplot(2,3,6);
% imhist(HH);
% title('Histograma de la imagen filtrada');
% xlabel('Intensidad de los pixeles');
% ylabel('Cantidad de pixeles');
