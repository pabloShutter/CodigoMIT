clear,clc,close all;

path(path, '../Library');

%% Configuraci�n del simulador
filt=1;   %'1' canal UAC, '0' canal ideal
noise=0;  %'0' sin ruido,'1' AWGN,'2' UWAN

%% Par�metros QPSK
Ps=30; % Potencia de la se�al
d=sqrt(2*Ps);

%% Par�metros OFDM
nSymb=50;         % N�mero de simbolos OFDM de la secuencia

fs=500e3;        % Frecuencia de muestreo
Tm=1/fs;         % Periodo de muestreo

N=4096*2;          % N�mero de puntos de la FFT
M=N*3/4;          % N�mero de muestras del prefijo c�clico
Nutil=784*2;      % Numero de portadoras �tiles

T=(N+M)*nSymb/fs;      % Duracion de la simulacion
t=0:Tm:T-Tm;     % eje de tiempos
Fsimb=1/T;       % Frec de s�mbolo
Df=fs/N;             % Separaci�n entre portadoras
f=-fs/2:Df:fs/2-Df;  % eje de frecuencia OFDM

fc=80e3;    % Frecuencia portadora modulaci�n paso banda (DMT)
% 
% portTx=0;   % Portadora transmitida: puede tomar valores entre 1 y Nutil
%             % para tx una �nica portadora o '0' para tx todas las portadoras

%% Parametros del Canal
if noise==1
    N0_2=7e-14;   % DEP de ruido AWGN
elseif noise==2
    N0_2=3.5e-14;   % DEP de ruido UWAN
end

% Par�metros del modelo UAC
fmin=20e3;
fmax=140e3;
d0=100;   % Separaci�n Tx-Rx (m)
wTx=14;   % profundidad del transmisor hasta el fondo
wRx=14;   % profundidad del receptor hasta el fondo
w=20;     % Profundidad del agua desde la superficie
pmax=50;  % N�mero m�ximo de rayos significativos considerados
dib=0;    % Representaci�n respuesta del canal: '0'->No, '1'->S�