% -----------------------------------------------------------
% freqIndex.m
%
% Devuelve el índice Matlab (semieje negativo a continuación
% del positivo) para una determinada frecuencia que se le pasa
% como parámetros. Si esa frecuencia concreta no está en el eje
% de frecuencia, se devuelve el índice de la siguiente frecuencia
% más cercana.
% -----------------------------------------------------------
% function [index,f]=freqIndex(freq,fs,Df)
% -----------------------------------------------------------
% Entrada/s:
%   freq: valor de la frecuencia de entrada
%   fs:   frecuencia de muestreo
%   Df:   paso en frecuencia
% Salidas:
%   index: índice de la frecuencia 'freq' o una frec cercana
%   f:     frec correspondiente a 'index' en el eje de frecuencias
% -----------------------------------------------------------

function [index,f]=freqIndex(freq,fs,Df)

%eje de frecuencia Matlab
tam=int32(fs/Df);
f_axis=zeros(1,tam);
f_axis(1:tam/2)=0:Df:fs/2-Df;
f_axis(tam/2+1:end)=-fs/2:Df:-Df;

index=find(f_axis>=abs(freq),1);
if freq<0
    index=(tam+1)-(index-1);
end
f=f_axis(index);
