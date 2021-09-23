%  function [a,bb]=Decisor2DPSK(q)
%  para la decisi�n de s�mbolos 2DPSK 
%  q    vector de muestras demoduladas en plano complejo
%
%  a     secuencia de s�mbolos decididos
%  bb    secuencia de bits decididos

function [a,bb]=DecisorDQPSK(q,A)
N=length(q);    % n�mero de s�mbolos

rho=abs(A/2+(A/2)*1j);

a=zeros(1,N);
bb=zeros(1,N);

a(1)=pol2cart(0, rho);
bb(1)=1;

for i=2:N    
    difFase=angle(q(i)*conj(q(i-1)));
    
    if difFase>deg2rad(-45) && dif<=deg2rad(45)        %d=(1,1)
        l1=dot(q(i),q(i-1));
        if l1>0
            bb(i)=1;
            a(i)=a(i-1);
        else
            bb(i)=0;
            a(i)=-a(i-1);
        end
    elseif difFase>deg2rad(45) && dif<=deg2rad(135)    %d=(0,1)
        
    elseif (difFase>deg2rad(135) && difFase<=deg2rad(180)) || (dif>=deg2rad(-180) && dif<=deg2rad(-135)) %d=(0,0)
        
    elseif difFase>deg2rad(-135) && dif<=deg2rad(-45)  %d=(1,0)
        
    end                    
end