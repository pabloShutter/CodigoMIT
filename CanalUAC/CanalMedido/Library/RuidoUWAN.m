function n=RuidoUWAN(t,N0_2,dib)

fm=1/(t(2)-t(1));
s=0.4;    %nivel de actividad de barcos entre 0 y 1
w=3.6;    %velocidad del viento (m/s)
f=[1 50:50:1e5]; %eje de frecuencia en el que se calcula de DEP del modelo UWAN

% 1. Creo el ruido AWG con Sw=0.2*10-9 W/Hz
potencia = N0_2 * fm;  %Potencia = N0 * fm/2 = N02 * fm
awgn = sqrt(potencia) * randn(1,length(t));

%2. Calculo DATO(f): la DEP de ruido ambiental
DEP_mod=UAC_LTI_Na(f,s,w); %d.e.p en (uPa)^2/Hz

if dib
    figure, semilogx(f,10*log10(DEP_mod)), grid on, xlabel('f (Hz)')
    ylabel('10log10(DATO(f)) dB re (1µPa)^2/Hz')
    title('d.e.p. del ruido acústico ambiental, s=0.4 y w=3.6 m/s')
end

%3. Calculo la respuesta al impulso, hn[n], del filtro
Hn=sqrt(DEP_mod*N0_2^-1); %resp en frec
h0=ifft(Hn);
energiah0=h0*h0';
hn=h0/sqrt(energiah0); %resp al impulso normalizada

if dib
    figure, subplot(2,1,1), semilogx(f,10*log10(abs(Hn))),ylabel('10log10(|Hn|) dB')
    xlabel('f(Hz)'), grid on, title('Respuesta en frecuencia del filtro')
    subplot(2,1,2), plot(abs(hn)), grid on, title('Respuesta al impulso del filtro')
    ylabel('|hn|')
end

%4. Filtro el ruido AWG con esta hn[n]
n=real(filter(hn,1,awgn));
% n=n(round(length(hn)/2):end); 
% n=n(1:end-floor(length(hn)/2));

if dib
    figure
    DEP_n=EstimaDEP(n,length(f),fm,'b');
    
    mitad=round(length(DEP_n)/2);
    w=pi*[0:1/mitad:1];
    freq=w*fm/(2*pi);

    figure, subplot(2,1,1), plot(t,real(n)), grid on
    xlabel('t(s)'),ylabel('n(t)'), title('Ruido gaussiano y coloreado')
    subplot(2,1,2), plot(freq/1e3,10*log10(DEP_n(1:mitad+1))),ylabel('10log10(DEP_n) dB'), xlabel('f (KHz)'),grid on
    title('d.e.p. del ruido generado')
end