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
    vent4=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el símbolo OFDM
    x=xx(vent4);
    
    xo=x((1+M):(N+M));
    Qx=RxOFDM(xo,FEQ,Nutil,iFc);
    
    ss=TxOFDM(Qx,N,0,Nutil,iFc,iFcH);
    vent3=(kSymb-1)*N+(1:N); % Ventana para el vector de símbolos OFDM
    yy(vent3)=ss;
end