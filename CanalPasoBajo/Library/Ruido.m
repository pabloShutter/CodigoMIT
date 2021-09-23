%function y = Ruido(t,N02)
%  t              Vector con instantes de tiempo
%  N02            Densidad Espectral de Potencia del Ruido en W/Hz
%
%  Crea un señal de ruido blanco gaussiano y ergódico de la densidad espectral de potencia especificada y media nula      
% 
function y = Ruido(t,N02)
  
 fm = 1/(t(2)-t(1));
 potencia = N02 * fm;  %Potencia = N0 * fm/2 = N02 * fm
 
 y = sqrt(potencia) * randn(1,length(t));