%% Configuracion
noise=1;   % Ruido canal -> '0': sin ruido, '1': AWGN, '2': UWAN
filt=0;    % Canal con resp al impulso h

%% Parametros OFDM
nSymb=300;       % N�mero de simbolos OFDM de la secuencia

fs=500e3;        % Frecuencia de muestreo
Tm=1/fs;         % Periodo de muestreo

a=1;
N=a*1024;          % N�mero de puntos de la FFT
M=N*1/4;          % N�mero de muestras del prefijo c�clico
Nutil=a*196;       % Numero de portadoras �tiles

T=(N+M)*nSymb/fs;      % Duracion de la simulacion
t=0:Tm:T-Tm; % eje de tiempos
Fsimb=1/T;       % Frec de s�mbolo
Df=fs/N;             % Separaci�n entre portadoras
f=-fs/2:Df:fs/2-Df;  % eje de frecuencia OFDM

fc=80e3;    % Frecuencia portadora modulaci�n paso banda (DMT)

%QPSK
Ps=2;%5000;           % pot se�al tx
d=sqrt(2*Ps);   % dist m�nima entre s�mbolos de la constelaci�n QPSK

%% Parametros del Canal
N0_2=2e-10;
% SNRdB=8;
fcut=120e3;            % Frecuencia de corte del filtro FIR paso bajo (canal)
h=fir1(32,2*fcut/fs);    % Resp al impulso del filtro FIR
H=freqz(h,1,N,'whole');  %Calcula Resp. frec. canal
[m,sincro]=max(h);   % calculamos la muestra optima para luego sincronizar
ret=length(h)-1;

%representaci�n gr�fica respuesta canal
% Th=length(h)/fs; %duraci�n de la respuesta impulsiva del canal
% th=0:Tm:Th-Tm;       %eje de tiempo respuesta canal
% figure, subplot(2,1,1), stem(th*1e6,abs(h)), grid on, xlabel('t (\mus)'),title('|h(t)|')
% subplot(2,1,2), plot(f/1e3, abs(fftshift(H))), grid on, xlabel('f (KHz)'),title('|H(f)|')