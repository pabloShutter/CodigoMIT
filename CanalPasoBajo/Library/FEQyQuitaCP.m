% -----------------------------------------------------------
% quitaCP.m
%
% Recibe una se�al x compuesta por una secuencia de s�mbolos
% OFDM con sus extensiones c�clicas (CP) y devuelve esa misma
% se�al habiendo eliminado el CP y con el FEQ aplicado
% -----------------------------------------------------------
% function y = FEQyQuitaCP(x,N,M,Ga,nSymb)
% -----------------------------------------------------------
% Entradas: 
%   x: se�al OFDM (longitud: (N+M)*nSymb)
%   N: Numero de puntos de la FFT
%   M: Numero de muestras del prefijo c�clico
%   Nutil: Numero de portadoras utiles (simb. QPSK/simb OFDM)
%   Ga: FEQ
%   nSymb: Numero de simbolos OFDM de la secuencia
% -----------------------------------------------------------

function yy = FEQyQuitaCP(xx,N,M,Nutil,FEQ,nSymb,iFc,iFcH)
  
yy=zeros(1,N*nSymb); %Solo se�al sin ruido, compensada y sin CP

for kSymb=1:nSymb    
    ventOFDM=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el s�mbolo OFDM
    x=xx(ventOFDM);
    
    xo=x((1+M):(N+M)); % Se le quita la parte del CP a la se�al
    Qx=RxOFDM(xo,FEQ,Nutil,iFc);  % Se le aplica el FEQ en el dominio de f  
    
    ss=TxOFDM(Qx,N,0,Nutil,iFc,iFcH);  % Se pasa al dominio de t
    vent=(kSymb-1)*N+(1:N); % Ventana para el vector de s�mbolos OFDM
    yy(vent)=ss;
end