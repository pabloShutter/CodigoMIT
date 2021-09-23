% -----------------------------------------------------------
% Ruido AWGN
% Devuelve un vector temporal de ruido blanco y gaussiano con una DEP
% que se le pasa como parámetro
% -----------------------------------------------------------
% function y = Ruido(t,N02)
% -----------------------------------------------------------
% Entradas:
%   t       Vector con instantes de tiempo
%   N02     Densidad Espectral de Potencia del Ruido en W/Hz
% -----------------------------------------------------------
function y = Ruido(t,N02)
  
 fs = 1/(t(2)-t(1));   % Frecuencia de muestreo
 Pn = N02 * fs;  % Potencia de ruido
 
 y = sqrt(Pn) * randn(1,length(t));