load('C:\Users\pablo\OneDrive - Universidad de M�laga\TFM_MIT\Simulador\CANAL_UAC\CanalUACmedido\LTVch_medidos\2017-06-02 09.47.44_ME_D50_S12.5__PB__SMultitono_fo80k_B96k_N1609_CAL.mat')
% load('C:\Users\pablo\OneDrive - Universidad de M�laga\TFM\Simulador\CANAL_UAC\CanalUACmedido\LTVch_medidos\RespuestasAlineadasRayoPpal\2017-06-02 09.47.44_ME_D50_S12.5__SMultitono_fo80k_B96k_N1609__ALIN.mat')

%% Configuraci�n del simulador
dib=0;
ch=1;   %'1' canal UAC, '0' canal ideal
tipoRuido=1;  %'0' sin ruido,'1' AWGN,'2' UWAN

%% Par�metros QPSK
Ps=2; % Potencia de la se�al
d=sqrt(2*Ps);

%% Par�metros OFDM
nSymb=300;         % N�mero de simbolos OFDM de la secuencia

fs=500e3;        % Frecuencia de muestreo
Tm=1/fs;         % Periodo de muestreo

N=4096*2;          % N�mero de puntos de la FFT
M=N*3/4;          % N�mero de muestras del prefijo c�clico
Nutil=600*2;      % Numero de portadoras �tiles

Df=fs/N;             % Separaci�n entre portadoras

fc=80e3;    % Frecuencia portadora modulaci�n paso banda (DMT)
%% Parametros del Canal
if tipoRuido==1
    N0_2=5e-14;   % DEP de ruido AWGN
elseif tipoRuido==2
    N0_2=5e-15;   % DEP de ruido UWAN
end
% N0_2=1e-10;   % DEP de ruido

% Respuesta variante con el tiempo (LTV) del canal UAC medido
canal=1;  % Receptor seleccionado (se utilizaron dos receptores)
h=TimeResponse.h{canal};   % Resp imp LTV
Delay=TimeResponse.delay(TimeResponse.delay<=10e-3);
t=TimeResponse.t;

Tm=Delay(2)-Delay(1); %Periodo de muestreo
fm=1/Tm;              %Frecuencia de muestreo

%Representaci�n resp imp LTV mapa de colores
if dib==1
    figure; fig=gcf;
    set(fig,'units','normalized','outerposition',[0 0 1 1],'name','Respuesta al impulso variante');
    color={'r','g','b','c','k','y'};
    colormap('jet');
    imagesc(Delay*1e3, t,20*log10(abs(h(:,Delay<=10e-3))) - 20*log10(max(max(abs(h)))),[-40,0])

    colorbar; set(gca, 'YDir','normal')
    xlabel(' delay (ms) '); ylabel(' time (sec) ');
    title( [' Respuesta impulsiva variante del canal ', num2str(canal)] );
end