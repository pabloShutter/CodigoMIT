%  ---------------------------------------------------------------------------
%  function [a,bb]=DecisorDQPSK(q)
%  ---------------------------------------------------------------------------
%  QPSK diferencial para la decisi�n de s�mbolos QPSK recibidos
%
%  NOTA: es importante que en el transmisor se env�e un primer s�mbolo QPSK
%  de referencia de valor q=1+1j
%  ---------------------------------------------------------------------------
%  q    vector de muestras demoduladas en plano complejo
%
%  a       secuencia de s�mbolos decididos
%  bb      secuencia de bits decididos
%  ---------------------------------------------------------------------------

function [a,bb]=DecisorDQPSK(q,d)
tam=length(q);    % n�mero de s�mbolos

%Primer s�mbolo QPSK de referencia (1,1)
a(1)=(d/2)+(d/2)*1j; % s�mbolo estimado
b(1,:)=[1,1];        % bits estimados en matriz Nx2

% creamos las secuencias de salida del decisor
for i=2:tam
    
   dif=angle(q(i-1)*conj(q(i))); %diferencia de fase entre s�mbolo actual y anterior
   
   if dif<=deg2rad(45) && dif>deg2rad(-45) %mismo cuadrante que el s�mbolo anterior
       if a(i-1)==((d/2)+(d/2)*1j) % simbolo anterior es (1,1)
           a(i)=(d/2)+(d/2)*1j;
           b(i,:)=[1,1];
       elseif a(i-1)==(-(d/2)+(d/2)*1j) % simbolo anterior es (0,1)
           a(i)=-(d/2)+(d/2)*1j;
           b(i,:)=[0,1];
       elseif a(i-1)==(-(d/2)-(d/2)*1j) % simbolo anterior es (0,0)
           a(i)=-(d/2)-(d/2)*1j;
           b(i,:)=[0,0];
       elseif a(i-1)==((d/2)-(d/2)*1j) % simbolo anterior es (1,0)
           a(i)=(d/2)-(d/2)*1j;
           b(i,:)=[1,0];
       end
       
   elseif dif<=deg2rad(-45) && dif>deg2rad(-135)
       if a(i-1)==((d/2)+(d/2)*1j) % simbolo anterior es (1,1)
           a(i)=-(d/2)+(d/2)*1j;
           b(i,:)=[0,1];
       elseif a(i-1)==(-(d/2)+(d/2)*1j) % simbolo anterior es (0,1)
           a(i)=-(d/2)-(d/2)*1j;
           b(i,:)=[0,0];
       elseif a(i-1)==(-(d/2)-(d/2)*1j) % simbolo anterior es (0,0)
           a(i)=(d/2)-(d/2)*1j;
           b(i,:)=[1,0];
       elseif a(i-1)==((d/2)-(d/2)*1j) % simbolo anterior es (1,0)
           a(i)=(d/2)+(d/2)*1j;
           b(i,:)=[1,1];
       end
       
   elseif (dif<=deg2rad(-135) && dif>=deg2rad(-180)) || (dif<=deg2rad(180) && dif>deg2rad(135))
       if a(i-1)==((d/2)+(d/2)*1j) % simbolo anterior es (1,1)
           a(i)=-(d/2)-(d/2)*1j;
           b(i,:)=[0,0];
       elseif a(i-1)==(-(d/2)+(d/2)*1j) % simbolo anterior es (0,1)
           a(i)=(d/2)-(d/2)*1j;
           b(i,:)=[1,0];
       elseif a(i-1)==(-(d/2)-(d/2)*1j) % simbolo anterior es (0,0)
           a(i)=(d/2)+(d/2)*1j;
           b(i,:)=[1,1];
       elseif a(i-1)==((d/2)-(d/2)*1j) % simbolo anterior es (1,0)
           a(i)=-(d/2)+(d/2)*1j;
           b(i,:)=[0,1];
       end
   elseif dif<=deg2rad(135) && dif>deg2rad(45)
       if a(i-1)==((d/2)+(d/2)*1j) % simbolo anterior es (1,1)
           a(i)=(d/2)-(d/2)*1j;
           b(i,:)=[1,0];
       elseif a(i-1)==(-(d/2)+(d/2)*1j) % simbolo anterior es (0,1)
           a(i)=(d/2)+(d/2)*1j;
           b(i,:)=[1,1];
       elseif a(i-1)==(-(d/2)-(d/2)*1j) % simbolo anterior es (0,0)
           a(i)=-(d/2)+(d/2)*1j;
           b(i,:)=[0,1];
       elseif a(i-1)==((d/2)-(d/2)*1j) % simbolo anterior es (1,0)
           a(i)=-(d/2)-(d/2)*1j;
           b(i,:)=[0,0];
       end
   end
end

bb=reshape(b,1,2*tam);   % convertimos matriz Nx2 en vector