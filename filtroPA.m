% FILTRO PASA ALTOS
% =================

% 1) IMPORTAMOS LA IMAGEN
G=imread('lenna.jpg');
G=rgb2gray(G);      % Sólo si la imagen es ".jpg"
GD=double(G);

% 2) TRANSF. FOURIER DE LA IMAGEN
GD=fft2(GD);        % Transformo la imagen
tfa=fftshift(GD);   % Normalizo
tfaabs=abs(tfa);    % Extraigo valores absolutos

% 3) VENTANA
[d1,d2]=size(GD);
ventana=ones(size(GD));
a=18;               % Base del obstáculo
b=4;                % Altura del obstáculo
ventana((d1/2)-a:(d1/2)+a,(d2/2)-b:(d2/2)+b)=0;

% 4) ACCIÓN DEL FILTRO
tfafilt=tfa.*ventana;
tfafiltabs=abs(tfafilt);

% 5) ANTITRANSF. FOURIER DE LA IMAGEN FILTRADA
afilt=abs(ifft2(tfafilt));
maxafilt=max(max(afilt));
H=afilt/maxafilt;

% 6) GRAFICA
figure;
subplot(1,2,1);
imshow(G); title('Imagen original');
subplot(1,2,2);
imshow(H); title('Imagen filtrada');
