%  ---------------------------------------------------------------------------
%  function [a,bb]=DecisorDBPSK(q)
%  ---------------------------------------------------------------------------
%  QPSK diferencial para la decisi�n de s�mbolos BPSK recibidos
%
%  NOTA: es importante que en el transmisor se env�e un primer s�mbolo BPSK
%  de referencia de valor q=1
%  ---------------------------------------------------------------------------
%  q    vector de muestras demoduladas en plano complejo
%
%  a       secuencia de s�mbolos decididos
%  bb      secuencia de bits decididos
%  ---------------------------------------------------------------------------

function [a,b]=DecisorDBPSK(q,d,BPSKtx)
N=length(q);    % n�mero de s�mbolos

a(1)=d/2+1j*d/2; b(1)=1;  
% creamos las secuencias de salida del decisor
for i=2:N
    
    dif=angle(q(i-1)*conj(q(i))); %diferencia de fase entre s�mb rx y s�mb de ref
    
    if dif<=deg2rad(90) && dif>deg2rad(-90) %est� en la misma zona de decisi�n que el s�mb ant
        a(i)=a(i-1);    % s�mbolos estimados
        if a(i-1)==((d/2)+(d/2)*1j)            
            b(i)=1;      % bits estimados
        elseif a(i-1)==(-(d/2)-(d/2)*1j)
            b(i)=0;
        end
    else
        a(i)=-a(i-1);    % s�mbolos estimados
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
           fprintf('Error: s�mb decidido a(%d) = %f + j%f  y  s�mb transmit A(%d) = %f + j%f\n',...
           i,real(a(i)),imag(a(i)),i,real(BPSK2(i)),imag(BPSK2(i))) 
       end
    end
end