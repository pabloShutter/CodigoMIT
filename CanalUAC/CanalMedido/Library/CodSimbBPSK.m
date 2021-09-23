% -----------------------------------------------------------
% Codificador BPSK
% Codificaci�n de bits como s�mbolos BPSK
% -----------------------------------------------------------
% function y=CodSimbQPSK(x)
% -----------------------------------------------------------
% x:  Secuencia de bits {0,1}
% y:  Secuencia de s�mbolos BPSK {+-1}
% -----------------------------------------------------------

function y=CodSimbBPSK(x,d)
tam=length(x);
y=[]; %Creamos vac�a la secuencia de s�mbolos BPSK

for i=1:tam
    if x(i)==1
        y(i)=d/2+1j*d/2;
    elseif x(i)==0
        y(i)=-d/2-1j*d/2;
    end
end