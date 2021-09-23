% -----------------------------------------------------------
% Codificador QPSK
% Codificación de bits como símbolos QPSK
% -----------------------------------------------------------
% function y=CodSimbQPSK(x)
% -----------------------------------------------------------
% x:  Secuencia de bits {0,1}
% y:  Secuencia de símbolos complejos {+-1,+-j}
% -----------------------------------------------------------

function y=CodSimbDQPSK(x,A)
N=length(x);
xx=reshape(x,N/2,2); % Organizamos los bits en 2 columnas

y=[]; %Creamos vacía la secuencia de símbolos QPSK (2 bits/simb)

rho=abs(A/2+(A/2)*1j); %Módulo del símbolo DQPSK transmitido

d=zeros(N/2+1,2);

d(1,:)=[1,1]; %Estado inicial

for i=2:N/2+1
    d(i,1)=~xor(xx(i-1,1),d(i-1,1));
    d(i,2)=~xor(xx(i-1,2),d(i-1,2));
    
    if d(i,:)==[1,1]       %d=(1,1)
        [re,im]=pol2cart(0, rho);
        y(i-1)=re+im*1i;
    elseif d(i,:)==[0,1]   %d=(0,1)
        [re,im]=pol2cart(pi/2, rho);
        y(i-1)=re+im*1i;
    elseif d(i,:)==[0,0]   %d=(0,0)
        [re,im]=pol2cart(pi, rho);
        y(i-1)=re+im*1i;
    elseif d(i,:)==[1,0]   %d=(1,0)
        [re,im]=pol2cart(3*pi/2, rho);
        y(i-1)=re+im*1i;
    end
end