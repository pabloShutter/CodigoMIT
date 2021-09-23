% -----------------------------------------------------------
% Codificador 2DPSK (2PSK Diferencial)
% Codificación de bits como símbolos 2SPSK
% -----------------------------------------------------------
% function y=CodSimb2DPSK(x)
% -----------------------------------------------------------
% x:  Secuencia de bits {0,1}
% y:  Secuencia de símbolos complejos {+-1,+-j}
% -----------------------------------------------------------

function y=CodSimb2DPSK(x,A)
N=length(x);

b=zeros(1,N+1);
b(2:end)=x;

y=zeros(1,N+1); %Creamos vacía la secuencia de símbolos QPSK

d=zeros(1,N+1);    

fase=zeros(1,N+1); %fase transmitida
rho=abs(A/2+(A/2)*1j);

%Estado inicial
d(1)=1;
fase(1)=0;
y(1)=pol2cart(fase(1), rho);

for i=2:N+1
    d(i)=~xor(b(i),d(i-1));
    
    if(d(i)==1)
        fase(i)=0;
    else
        fase(i)=pi;
    end
    
    y(i)=pol2cart(fase(i), rho);
end

y=y(2:end); %Nos deshacemos del símbolo de referencia porque el receptor lo 
            %conoce de antemano