% -----------------------------------------------------------
% Receptor OFDM
% Elimina el prefijo ciclico
% Calcula el símbolo en el dominio de la frecuencia
% Reordena los símbolos (deshace la reordenación del transmisor)
% -----------------------------------------------------------
% function V=RxOFDM12(v,N,M,Nutil,dib)
% -----------------------------------------------------------
% v: simbolo de entrada (longitud: N+M)
% V: vector de simbolos de salida (longitud: Nutil)
% N: Numero de puntos de la FFT
% M: Numero de muestras del prefijo cíclico
% Nutil: Numero de portadoras útiles
% dib: 1: dibujar gráficas, 0: no dibujar
% -----------------------------------------------------------

function Q=RxOFDM(v,FEQ,Nutil,iFc)

Yo=fft(v);                  % FFT

Y=Yo.*FEQ; %Compensación del canal

% Deshago la interpolación en frecuencia para quedarme con los Nutil
% símbolos QPSK:
Q=zeros(1,Nutil);
Q(1:Nutil/2)=Y(iFc-Nutil/2:iFc-1);
Q(Nutil/2+1:Nutil)=Y(iFc+1:iFc+Nutil/2);