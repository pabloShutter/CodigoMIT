%Actualización 14/06/2017 para escritura tesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Optimización de tiempos de
            % Función para evaluar la H(f) del modelo UAC LTI aguas someras
            % J.F.Paris
            % 13-06-2013
        %por Adrián S.G.
        %20-05-2015
        %Corrección de errores en cálculo de ángulos y distancias
            %Sólo se consideraban los caminos con número impar de rebotes.
            %Incorporamos los caminos de número par de rebotes.
            %También está mal orientado el diagrama de radiación de los
            %transductores. 
            %Estaba como si los transductores miraran hacia arriba,
            %cuando deben mirar hacia abajo.
        %por Adrián S.G.
        %27-01-2016
        %Corrección para permitir distintas profundidades de Tx y Rx
            %Ahora se permite que transductores estén a distinta profundidad
        %por Adrián S.G.
        %01-02-2016
        %Corrección error coeficiente de reflexión
            %Al calcular el coeficiente de reflexión en el fondo se usaba la
            %ecuación de Stojanovic y Qarabaqi, la cual usa el ángulo de incidencia
            %sobre el receptor (grazing angle), y se estaba usando con el ángulo de incidencia
            %sobre el fondo (ángulo de incidencia sobre el fondo = pi/2 - grazing angle)
            %Ya está corregido.
        %por Adrián S.G.
        %02-02-2016
        %Corrección error coeficiente de reflexión
            %Realmente cuando aux<0, debe salir Gamma compleja, con módulo
            %igual a 1 pero con una fase.
            %Ya está corregido.
        %por Adrián S.G.
        %02-02-2016
        %Inclusión de posibilidad de elegir modelo de atenuación de Marsh & Schulkin
            %Nuevos parámetros de entrada: Salinidad y temperatura del agua
        %por Adrián S.G.
        %13/06/2017
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Entrada/s:
        %   f       frecuencia en Hz
        %   d0      distancia TX-RX en m
        %   w       profundidad del agua
        %   wTx     profundidad del transmisor hasta el fondo
        %   wRx     profundidad del receptor hasta el fondo
        %   aTX     apertura vertical del diag. polar del TX en rad
        %           aTX=[aTX1 aTX2]
        %   aRX     apertura vertical del diag. polar del RX en rad
        %           aRX=[aRX1 aRX2]
        %   cb      velocidad de propagación del sonido en el fondo (m/s)
        %   rhob    densidad volumétrica del fondo (kg/m3)
        %   modelo  modelo de atenuación dependiente de la frecuencia
        %           puede ser 'thorp' (Thorp, 100Hz-3kHz),
        %           'm&s' (Marsh & Schulkin, 3kHz-500kHz),
        %           o 'bi' thorp y m&s por tramos
        %   S       salinidad del agua en tanto por mil, para modelo 'm&s'
        %   T       tempoeratura del agua en ºC, para modelo 'm&s'
        % Salida/s:
        %   H(f)    respuesta en frecuencia, presión RX / presión a 1m en TX 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rutina principal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H=UAC_LTI_H(f,d0,wTx,wRx,w,aTX,aRX,cb,rhob,pmax,modelo,S,T)

% parámetros libres del modelo del modelo

r=1.5;         % factor de spreading
c=1500;        % velocidad de propagación del sonido en el agua (m/s)
rho=1000;      % densidad volumétrica del agua (kg/m3)
ate_ex=40;     % nivel atenuación relativo para aceptación de un rayo (dB)


lp=0:pmax;  % lista potencial de rayos
%Los rayos se entienden de la siguiente forma
    %p=0, LOS
    %p impar, primer rebote en la superficie.
        %Número total de rebotes: b = (p+1)/2
    %p par, primer rebote en el fondo
        %Número total de rebotes: b = p/2


th=theta(lp); %th(lp+1) para lp>0 es el ángulo de indicencia del rayo sobre las superficies, es decir entre el rayo y la normal a la superficie
              % El ángulo th en tal caso sería el medido donde se marca el símbolo & en el esquema:
                % _________________________________
                %      |   |    /&|&\    |   |
                %     (Tx) |   /  |  \   |  (Rx)
                %       \  |  /   |   \  |  /
                % _______\&|&/____|____\&|&/_______
                %
              %En el caso de th(lp+1) para lp=0, el ángulo al que nos referimos sería el grazing angle en el receptor
              %El ángulo th en tal caso sería el medido donde se marca el símbolo & en el esquema.
              %En el ejemplo del esquema este th sería positvio, en caso de wRx>wTx, este sería negativo
                % ___________________
                %      |        |
                %     (Tx)      |
                %         \     |
                %          \    |
                %           \   |
                %            \  |
                %            &\ |
                %      --------(Rx)
                % ___________________
                
gam=gamma(lp); %Coeficientes acumulado de reflexión de cada rayo
dist=d(lp); %Distancia recorrida por cada rayo

alpha=ate(f);

Coef=gam./sqrt(dist.^r);


for p=lp+1,
    
    hp(p,:)=Coef(p)./sqrt(alpha.^dist(p));
    
end
    
atenuaciones=20*log10(abs(hp));
ate_rel=repmat(max(atenuaciones,[],1),pmax+1,1)-atenuaciones;
    
% filtrado de la lista potencial de rayos
lista_p=zeros(size(atenuaciones));

for p=lp,
    if p==0
        if (-th(p+1)<=(aTX(2)-pi/2)) && (-th(p+1)>=(aTX(1)-pi/2)) && (th(p+1)<=(aRX(2)-pi/2)) && (th(p+1)>=(aRX(1)-pi/2))
            lista_p(p+1,ate_rel(p+1,:)<ate_ex)=1;
        end
    elseif (rem(p,4)==1)
        if (th(p+1)>=(pi-aTX(2))) && (th(p+1)<=(pi-aTX(1))) && (th(p+1)>=(pi-aRX(2))) && (th(p+1)<=(pi-aRX(1)))
            lista_p(p+1,ate_rel(p+1,:)<ate_ex)=1;
        end
    elseif (rem(p,4)==2)
        if (th(p+1)<=aTX(2)) && (th(p+1)>=aTX(1)) && (th(p+1)<=aRX(2)) && (th(p+1)>=aRX(1))
            lista_p(p+1,ate_rel(p+1,:)<ate_ex)=1;
        end
    elseif (rem(p,4)==3)
        if (th(p+1)>=(pi-aTX(2))) && (th(p+1)<=(pi-aTX(1))) && (th(p+1)<=aRX(2)) && (th(p+1)>=aRX(1))
            lista_p(p+1,ate_rel(p+1,:)<ate_ex)=1;
        end
    else
        if (th(p+1)<=aTX(2)) && (th(p+1)>=aTX(1)) && (th(p+1)>=(pi-aRX(2))) && (th(p+1)<=(pi-aRX(1)))
            lista_p(p+1,ate_rel(p+1,:)<ate_ex)=1;
        end
    end
end


% cálculo final de la respuesta en frecuencia
fd=dist'*f;
H=sum(hp.*lista_p.*exp(-2*pi*1i*fd/c),1);




    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subrutinas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=theta(p)
    
    for k=1:length(p),
        if p(k)==0;  %Rayo directo
            y(k)=atan((wTx-wRx)/d0);
        elseif rem(p(k),2)==1  %Rayos con primera reflexión en superficie
            b=(p(k)+1)/2;
            if rem(b,2)==1
                y(k)=atan(d0/((b+1)*w-wTx-wRx));
            else
                y(k)=atan(d0/(b*w-wTx+wRx));
            end
        else                   %Rayos con primera reflexión en fondo
            b=p(k)/2;
            if rem(b,2)==1
                y(k)=atan(d0/((b-1)*w+wTx+wRx));
            else
                y(k)=atan(d0/(b*w-wRx+wTx));
            end
        end
    end
end % theta        
        

function y=d(p)

    for k=1:length(p),
        if p(k)==0;
            y(k)=d0*sec(th(k));
        elseif rem(p(k),2)==1
            b=(p(k)+1)/2;
            if rem(b,2)==1
                y(k)=((b+1)*w-wTx-wRx)*sec(th(k));
            else
                y(k)=(b*w-wTx+wRx)*sec(th(k));
            end
        else
            b=p(k)/2;
            if rem(b,2)==1
                y(k)=((b-1)*w+wTx+wRx)*sec(th(k));
            else
                y(k)=(b*w-wRx+wTx)*sec(th(k));
            end
        end
	end
    
end % distancias

function y=gammab(t)

    aux=(c/cb)^2-sin(t)^2;

    y=(rhob/rho*cos(t)-sqrt(aux))/(rhob/rho*cos(t)+sqrt(aux));
    

end % reflexion fondo

function y=gamma(p)
    
    for k=1:length(p),
        if p(k)==0,
            y(k)=1;
        elseif rem(p(k),2)==1,
            b=(p(k)+1)/2;
            if rem(b,2)==1
                y(k)=(-1)^((b+1)/2)*gammab(th(k))^((b-1)/2);
            else
                y(k)=(-1)^(b/2)*gammab(th(k))^(b/2);
            end
        else
            b=p(k)/2;
            if rem(b,2)==1
                y(k)=(-1)^((b-1)/2)*gammab(th(k))^((b+1)/2);
            else
                y(k)=(-1)^(b/2)*gammab(th(k))^(b/2);
            end
        end
    end

end % reflexion acumulada

function y=ate(f)
    
    fkHz=f/1000;
    if strcmp(modelo,'thorp')
        a= 0.11*fkHz.^2./(1+fkHz.^2) +44*fkHz.^2./(4100+fkHz.^2) +2.75/10000*fkHz.^2 +0.003;
    else
        if strcmp(modelo, 'm&s')
            A=2.34e-6; %Constante
            B=3.38e-6; %Constante

            g=9.8; %Aceleración gravedad
            prof=w-(wTx+wRx)/2; %Profundidad
            P=rho*g*prof; %Presión hidrostática en N/m2
            P=P*1e-5; %Presión en kg/cm2

            fT=21.9*10^(6-1520/(T+273));

            a= 8.68e3 * ( (S*A*fT*fkHz.^2)./(fT^2+fkHz.^2) + B*fkHz.^2/fT ) * (1-6.54e-4*P);
        else
            if strcmp(modelo,'bi')
                a= 0.11*fkHz.^2./(1+fkHz.^2) +44*fkHz.^2./(4100+fkHz.^2) +2.75/10000*fkHz.^2 +0.003;
                A=2.34e-6; %Constante
                B=3.38e-6; %Constante

                g=9.8; %Aceleración gravedad
                prof=w-(wTx+wRx)/2; %Profundidad
                P=rho*g*prof; %Presión hidrostática en N/m2
                P=P*1e-5; %Presión en kg/cm2

                fT=21.9*10^(6-1520/(T+273));

                b= 8.68e3 * ( (S*A*fT*fkHz.^2)./(fT^2+fkHz.^2) + B*fkHz.^2/fT ) * (1-6.54e-4*P);

                a= a.*(fkHz<=3) + b.*(fkHz>3);
            end
        end
    end
    y=10.^(a/10000);

end % atenuacíón

%EscribirCache(H);

end
