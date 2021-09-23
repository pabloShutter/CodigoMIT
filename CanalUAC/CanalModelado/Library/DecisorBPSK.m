%  function [a,bb]=DecisorBPSK(q)
%  para la decisión de símbolos BPSK 
%  q    vector de muestras demoduladas en plano complejo
%
%  a     secuencia de símbolos decididos
%  b    secuencia de bits decididos

function [a,b]=DecisorBPSK(q,d)
N=length(q);    % número de símbolos

a(1)=d/2+1j*d/2; b(1)=1;  
% creamos las secuencias de salida del decisor
for i=2:N
    
    dif=angle(q(i)*conj(q(1))); %diferencia de fase entre símb rx y símb de ref
    
    if dif<=deg2rad(90) && dif>deg2rad(-90)
       a(i)=d/2+1j*d/2;    % símbolos estimados
       b(i)=1;      % bits estimados
    else
       a(i)=-d/2-1j*d/2;
       b(i)=0;
    end
%     if real(q(i))>=0
%        a(i)=1;
%        b(i)=1;
%     elseif real(q(i))<0
%         a(i)=-1;    % símbolos estimados
%         b(i)=0;     % bits estimados
%     end    
end