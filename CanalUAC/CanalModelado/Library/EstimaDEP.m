function y=EstimaDEP(x,N,fm,color)
% function y=EstimaDEP(x,N,fm,color)
% Estima y Dibuja la DEP de una señal en forma de vector de muestras
% [Se utiliza el metodo de promediado de periodogramas enventanados con 
% una ventana de Hanning y no solapados.]
%
% x: vector de señal
% N: numero de puntos de la FFT
% fm: frecuencia de muestreo (en Hz)
% color: string para el color y trazo de la curva, ejemplo de valores
% posibles 'b','r','k','g'

tam=length(x); %tamaño señal a analizar
if ~(tam>N*2),
    disp('ERROR:'),disp('el tamaño del vector de señal es pequeño para el número de muestras de la FFT solicitado')
else if ~(tam>N*5),
        disp('AVISO:'),disp('el tamaño del vector de señal es pequeño para el número de muestras de la FFT solicitado y la estimación de la DEP no será buena')
    else
        win = .5*(1 - cos(2*pi*(1:N)/(N+1)));	%ventana de Hanning
        y=zeros(1,N);  %vector DEP
        z=zeros(1,N);  %vector periodogramas
        
        for k=1:(floor(tam/N)-1),   %se repite para tantos periodogramas como quepan
            z=x(k*N:(k+1)*N-1);   %selección del segmento de la señal
            y=y+(abs(fft(z.*win))).^2/sum(win.^2);
            %enventanado y periodograma, y acumulación
        end
        y=y/floor(tam/N); %normalización para tener el promedio: estima de la DEP discreta
	y=y/fm;	%normalización para pasar a DEP de señal continua
        
        mitad=round(N/2);      %solo se representa el eje de frecuencia positivo
        w=pi*[0:1/mitad:1];
        freq=w*fm/(2*pi);
        
%         plot(freq/1e3,(y(1:mitad+1)),color)
%         xlabel('frecuencia (kHz)'),
%         ylabel('DEP (dB)'),
%         grid on, zoom on;
    end
end

% Programa realizado por Francisco Javier Cañete
% Universidad de Málaga - Departamento de Ingeniería de Comunicaciones
