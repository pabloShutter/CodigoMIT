%  function [a,bb]=Decisor2DPSK(q)
%  para la decisión de símbolos 2DPSK 
%  q    vector de muestras demoduladas en plano complejo
%
%  a     secuencia de símbolos decididos
%  bb    secuencia de bits decididos

function [a,bb]=Decisor2DPSK(q,A)
N=length(q);    % número de símbolos

rho=abs(A/2+(A/2)*1j);

a=zeros(1,N);
bb=zeros(1,N);

a(1)=pol2cart(0, rho);
bb(1)=1;

for i=2:N
%     difFase=angle(q(i)*conj(q(i-1)));
    l=dot(q(i),q(i-1));
    if l>0
        bb(i)=1;
        a(i)=a(i-1);
    else
        bb(i)=0;
        a(i)=-a(i-1);
    end
    
end