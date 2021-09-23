
clear,clc,close all;

%% Configuración del simulador
ch=1;   %'1' canal UAC, '0' canal ideal
noise=0;  %'0' sin ruido,'1' AWGN,'2' UWAN

%% Parámetros QPSK
Ps=2; % Potencia de la señal
d=sqrt(2*Ps);

%% Parámetros OFDM
nSymb=200;         % Número de simbolos OFDM de la secuencia

fs=500e3;        % Frecuencia de muestreo
Tm=1/fs;         % Periodo de muestreo

N=8380;          % Número de puntos de la FFT
M=N*3/4;          % Número de muestras del prefijo cíclico
Nutil=600*2;      % Numero de portadoras útiles
% 
% T=(N+M)*nSymb/fs;      % Duracion de la simulacion
% t=0:Tm:T-Tm;     % eje de tiempos
% Fsimb=1/T;       % Frec de símbolo
Df=fs/N;             % Separación entre portadoras
% f=-fs/2:Df:fs/2-Df;  % eje de frecuencia OFDM

fc=80e3;    % Frecuencia portadora modulación paso banda (DMT)
%% Parametros del Canal
if noise==1
    N0_2=1e-9;   % DEP de ruido AWGN
elseif noise==2
    N0_2=5e-10;   % DEP de ruido UWAN
end
% N0_2=5e-11;   % DEP de ruido

% Respuesta estática del canal UAC medido
if ch==1
    [h,t,H,f]=respuestaUACestatica;
    ret=length(h)-1;
end
% 
% figure('Name','Respuesta UAC estática medida'), subplot(2,1,1), plot(abs(h)), grid on
% xlabel('Muestras'), ylabel('|h(m)|'), title('Respuesta al impulso realista')
% subplot(2,1,2), plot(f/1e3, abs(fftshift(H))), grid on
% xlabel('f (KHz)'), ylabel('|H(f)|'), title('Respuesta en frecuencia realista')