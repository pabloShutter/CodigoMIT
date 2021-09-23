function [h,t,H,f]=respuestaUACestatica

% Cargado de la respuesta impulsiva LTV medida
load('C:\Users\pablo\OneDrive - Universidad de Málaga\TFM_MIT\Simulador\CANAL_UAC\CanalUACmedido\LTVch_medidos\2017-06-02 09.47.44_ME_D50_S12.5__PB__SMultitono_fo80k_B96k_N1609_CAL.mat','TimeResponse')

h=TimeResponse.h{1}(1,:)*60;  % Selección de una respuesta impulsiva estática en un determinado instante de tiempo
t=TimeResponse.delay;      % Eje de tiempo de la respuesta impulsiva estática medida

fs=1/(t(2)-t(1)); %frec muestreo

H=fft(h);           % resp en frec
Df=fs/length(H);    % resolución frec (Hz)
f=-fs/2:Df:fs/2-Df; % eje de frecuencia (Hz)

figure, subplot(2,1,1), plot(t*1e3, abs(h)), grid on, xlabel('t (ms)'), ylabel('|h(t)|')
subplot(2,1,2), plot(f/1e3, abs(fftshift(H))), grid on, xlabel('f (KHz)'), ylabel('|H(f)|')