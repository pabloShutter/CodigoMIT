%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Función para evaluar la Na(f) del modelo UAC LTI
% J.F.Paris
% 4-11-2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Entrada/s:
%   f       frecuencia en Hz
%   s       nivel de actividad de barcos entre 0 y 1
%   w       velocidad del viento (m/s)
% Salida/s:
%   Na(f)    d.e.p en (uPa)^2/Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rutina principal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N=UAC_LTI_Na(f,s,w)

f=f/1000;

Lt=17-30*log10(f);
Ls=40+20*(s-0.5)+26*log10(f)-60*log10(f+0.03);
Lw=50+7.5*w^(1/2)+20*log10(f)-40*log10(f+0.4);
Lth=-15+20*log10(f);

N=10.^(Lt/10)+10.^(Ls/10)+10.^(Lw/10)+10.^(Lth/10);
