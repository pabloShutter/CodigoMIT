% -----------------------------------------------------------
% quitaCP.m
%
% Recibe una señal x compuesta por una secuencia de símbolos
% OFDM con sus extensiones cíclicas (CP) y devuelve esa misma
% señal habiendo eliminado el CP y con el FEQ aplicado
% -----------------------------------------------------------
% function y = FEQyQuitaCP(x,N,M,Ga,nSymb)
% -----------------------------------------------------------
% Entradas: 
%   x: señal OFDM (longitud: (N+M)*nSymb)
%   N: Numero de puntos de la FFT
%   M: Numero de muestras del prefijo cíclico
%   Nutil: Numero de portadoras utiles (simb. QPSK/simb OFDM)
%   Ga: FEQ
%   nSymb: Numero de simbolos OFDM de la secuencia
% -----------------------------------------------------------

function yy = FEQyQuitaCP(xx,N,M,Nutil,FEQ,nSymb,iFc,iFcH)
  
yy=zeros(1,N*nSymb); %Solo señal sin ruido, compensada y sin CP

for kSymb=1:nSymb    
    ventOFDM=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el símbolo OFDM
    x=xx(ventOFDM);
    
    xo=x((1+M):(N+M)); % Se le quita la parte del CP a la señal
    Qx=RxOFDM(xo,FEQ,Nutil,iFc);  % Se le aplica el FEQ en el dominio de f  
    
    ss=TxOFDM(Qx,N,0,Nutil,iFc,iFcH);  % Se pasa al dominio de t
    vent=(kSymb-1)*N+(1:N); % Ventana para el vector de símbolos OFDM
    yy(vent)=ss;
end