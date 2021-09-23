% -----------------------------------------------------------
% RuidoUWAN.m
%
% Devuelve un vector temporal de ruido UWAN con un DEP coloreada
% que se calcula a partir de una DEP plana N0_2 mediante un filtrado
% con la forma de la DEP del modelo UAC_LTI_Na
% -----------------------------------------------------------
% function n = RuidoUWAN(t,N0_2,dib)
% -----------------------------------------------------------
% Entradas: 
%   t:    vector temporal
%   N0_2: DEP del ruido AWGN en W/Hz
%   dib:  representar gráficamente los resultados '1' o '0'
% -----------------------------------------------------------

function n=RuidoUWAN(t,N0_2,dib)

fs=1/(t(2)-t(1));
s=0.4;    %nivel de actividad de barcos entre 0 y 1
w=3.6;    %velocidad del viento (m/s)
f=1:5:250e3;

% 1. Creo el ruido AWG con Sw=0.2*10-9 W/Hz
Pn = N0_2 * fs;
awgn = sqrt(Pn) * randn(1,length(t));

%2. Calculo DEPm(f): la DEP de ruido ambiental del modelo
DEPm=UAC_LTI_Na(f,s,w); %d.e.p en (uPa)^2/Hz

if dib
    figure, semilogx(f,10*log10(DEPm)), grid on, xlabel('f (Hz)')
    ylabel('10log10(DATO(f)) dB re (1µPa)^2/Hz')
    title('d.e.p. del ruido acústico ambiental, s=0.4 y w=3.6 m/s')
end

%3. Calculo la respuesta al impulso, hn[n], del filtro
Hn=sqrt(DEPm*N0_2^(-1)); %resp en frec
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
n=filter(hn,1,awgn);
% n=conv(hn,awgn);

if dib
    DEP_n=EstimaDEP(n,1024,fs,'b');

    Df=fs/(2*N); freq=0:Df:(fs/2)-Df;
    figure, subplot(2,1,1), plot(real(n)), grid on
    xlabel('n'),ylabel('n(n)'), title('Ruido gaussiano y coloreado')
    subplot(2,1,2), semilogx(freq/1e3,10*log10(DEP_n)),ylabel('10log10(DEP_n) dB'), grid on
    title('d.e.p. del ruido generado')
end