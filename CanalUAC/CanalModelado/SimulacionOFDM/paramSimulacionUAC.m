clear,clc,close all;
path(path, '../Library');

%% Parámetros QPSK
Ps=2; % Potencia de la señal
d=sqrt(2*Ps);

%% Parámetros OFDM
nSymb=300;         % Número de simbolos OFDM de la secuencia

fs=500e3;        % Frecuencia de muestreo
Tm=1/fs;         % Periodo de muestreo

N=4096;          % Número de puntos de la FFT
M=N*3/8;          % Número de muestras del prefijo cíclico
Nutil=784;      % Numero de portadoras útiles

T=(N+M)*nSymb/fs;      % Duracion de la simulacion
t=0:Tm:T-Tm;     % eje de tiempos
Fsimb=1/T;       % Frec de símbolo
Df=fs/N;             % Separación entre portadoras
f=-fs/2:Df:fs/2-Df;  % eje de frecuencia OFDM

fc=80e3;    % Frecuencia portadora modulación paso banda (DMT)
% 
% portTx=0;   % Portadora transmitida: puede tomar valores entre 1 y Nutil
%             % para tx una única portadora o '0' para tx todas las portadoras

%% Parametros del Canal

% Parámetros del modelo UAC
fmin=20e3;
fmax=140e3;
d0=100;  % Separación Tx-Rx (m)
wTx=14;   % Profundidad Tx
wRx=14;   % Profundidad Rx
w=20;    % Profundidad fondo marino (desde superficie)
pmax=50; % Número máximo de rayos significativos considerados
dib=0;   % Representación respuesta del canal: '0'->No, '1'->Sí