%  function [a,bb]=DecisorQPSK(q)
%  para la decisión de símbolos QPSK 
%  q    vector de muestras demoduladas en plano complejo
%
%  a     secuencia de símbolos decididos
%  bb    secuencia de bits decididos

function [a,bb]=DecisorQPSK(q,d)
N=length(q);    % número de símbolos

% creamos las secuencias de salida del decisor
for i=1:N;
    if (real(q(i))>=0 && imag(q(i))<0),
        a(i)=(d/2)-(d/2)*1j;      % símbolo estimado
        b(i,:)=[1,0];   % bits estimados en matriz Nx2
    elseif (real(q(i))<0 && imag(q(i))<0)
        a(i)=-(d/2)-(d/2)*1j;
        b(i,:)=[0,0];
    elseif (real(q(i))<0 && imag(q(i))>0),
        a(i)=-(d/2)+(d/2)*1j;
        b(i,:)=[0,1];
    else
        a(i)=(d/2)+(d/2)*1j;
        b(i,:)=[1,1];
    end
end
bb=reshape(b,1,2*N);   % convertimos matriz Nx2 en vector
