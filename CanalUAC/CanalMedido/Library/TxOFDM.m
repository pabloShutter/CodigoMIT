% -----------------------------------------------------------
% Transmisor OFDM
%
% Coloca simbolos en portadoras utiles ordenadamente
% desde la frecuencia más negativa hasta la más positiva
% En ambos extremos deja a cero las frecuencias "de guarda"
% Deja también a cero a la portadora cero (DC)
% Calcula el símbolo en el dominio del tiempo mediante IFFT
% e inserción del prefijo cíclico
% -----------------------------------------------------------
% function s=TxOFDM(A,N,M,dib)
% -----------------------------------------------------------
% s: simbolo de salida (longitud: N+M)
% A: vector de símbolos de entrada (longitud: Nutil)
% N: Numero de puntos de la FFT
% M: Numero de muestras del prefijo cíclico
% dib: 1: dibujar gráficas, 0: no dibujar
% -----------------------------------------------------------

function s=TxOFDM(A,N,M,Nutil,iFc,iFcH)

A_conj=conj(A);

% 1. Coloca los símbolos en portadoras útiles (simetría hermítica)
Ao=zeros(1,N);
Ao(iFc-Nutil/2:iFc-1)=A(1:Nutil/2);
Ao(iFc+1:iFc+Nutil/2)=A(Nutil/2+1:Nutil);

Ao(iFcH+Nutil/2:-1:iFcH+1)=A_conj(1:Nutil/2);
Ao(iFcH-1:-1:iFcH-Nutil/2)=A_conj(Nutil/2+1:Nutil);

% 2. Señal en dominio del tiempo e inserción del CP
s0=ifft(Ao);            % IFFT
s=[s0((N-M+1):N) s0];   % Inserción del CP