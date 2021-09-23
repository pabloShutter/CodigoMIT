function Pxx = EstimaDEP_Pablo(x,N,fm)

tam=length(x); %tamaño señal a analizar
if ~(tam>N*2),
    disp('ERROR:'),disp('el tamaño del vector de señal es pequeño para el número de muestras de la FFT solicitado')
else if ~(tam>N*5),
        disp('AVISO:'),disp('el tamaño del vector de señal es pequeño para el número de muestras de la FFT solicitado y la estimación de la DEP no será buena')
    else
        % Creación de ventana de Hanning y pre-cálculo de la escala
        win = 0.5*(1-cos(2*pi*(1:N)/(N-1))); % win=hann(N);
        s_win = norm(win,2)^2;

        % Inicialización de variable Pxx para almacenar todos los segmentos
        Pxx = zeros(N,2);
        for k=1:(floor(tam/N)-1)
            x1 = x(1+k*N:(k+1)*N); % Se toma el segmento k de la señal
            xw = x1.*win;          % Se multiplica por la ventana de Hanning            
            Pxx(:,k+1) = abs(fft(xw)).^2./s_win; % Se ajusta la magnitud y la escala
        end
        
        Pxx = mean(Pxx,2); % Promedio de todos los segmentos
        Pxx = Pxx./(fm);  % Se normaliza de acuerdo a la frec de muestreo
        
        Pxx = Pxx(1:length(Pxx)/2+1);
        % Preserve the energy (we just chopped-off half of spectrum), but do not 
        % modify first frequency bin as it is unique (mean value of a signal)
        % Pxx(2:end-1) = 2*Pxx(2:end-1);

        mitad=round(N/2);      %solo se representa el eje de frecuencia positivo
        w=pi*(0:1/mitad:1);
        freq=w*fm/(2*pi);

        plot(freq/1e3,(Pxx(1:mitad+1)))
        xlabel('frecuencia (kHz)'),
        ylabel('DEP (dB)'),
        grid on
    end
end