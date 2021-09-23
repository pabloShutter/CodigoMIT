% -----------------------------------------------------------
% Codificador BPSK
% Codificación de bits como símbolos BPSK
% -----------------------------------------------------------
% function y=CodSimbQPSK(x)
% -----------------------------------------------------------
% x:  Secuencia de bits {0,1}
% y:  Secuencia de símbolos BPSK {+-1}
% -----------------------------------------------------------

function y=CodSimbBPSK(x,d)
tam=length(x);
y=[]; %Creamos vacía la secuencia de símbolos BPSK

for i=1:tam
    if x(i)==1
        y(i)=d/2+1j*d/2;
    elseif x(i)==0
        y(i)=-d/2-1j*d/2;
    end
end