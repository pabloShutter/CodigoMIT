%  ---------------------------------------------------------------------------
%  function [a,bb]=DecisorDBPSK(q)
%  ---------------------------------------------------------------------------
%  QPSK diferencial para la decisión de símbolos BPSK recibidos
%
%  NOTA: es importante que en el transmisor se envíe un primer símbolo BPSK
%  de referencia de valor q=1
%  ---------------------------------------------------------------------------
%  q    vector de muestras demoduladas en plano complejo
%
%  a       secuencia de símbolos decididos
%  bb      secuencia de bits decididos
%  ---------------------------------------------------------------------------

function [a,b]=DecisorDBPSK(q,d,BPSKtx)
N=length(q);    % número de símbolos

a(1)=d/2+1j*d/2; b(1)=1;  
% creamos las secuencias de salida del decisor
for i=2:N
    
    dif=angle(q(i-1)*conj(q(i))); %diferencia de fase entre símb rx y símb de ref
    
    if dif<=deg2rad(90) && dif>deg2rad(-90) %está en la misma zona de decisión que el símb ant
        a(i)=a(i-1);    % símbolos estimados
        if a(i-1)==((d/2)+(d/2)*1j)            
            b(i)=1;      % bits estimados
        elseif a(i-1)==(-(d/2)-(d/2)*1j)
            b(i)=0;
        end
    else
        a(i)=-a(i-1);    % símbolos estimados
        if a(i-1)==((d/2)+(d/2)*1j)            
            b(i)=0;      % bits estimados
        elseif a(i-1)==(-(d/2)-(d/2)*1j)
            b(i)=1;
        end
    end   
end

if length(BPSKtx)>1200
    BPSK2=BPSKtx(1201:2400).';
    for i=1:length(a)
       if a(i)~=BPSK2(i)
           fprintf('Error: símb decidido a(%d) = %f + j%f  y  símb transmit A(%d) = %f + j%f\n',...
           i,real(a(i)),imag(a(i)),i,real(BPSK2(i)),imag(BPSK2(i))) 
       end
    end
end