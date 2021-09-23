% -----------------------------------------------------------
% Canal LTV
%
% Toma como parámetro de entrada un símbolo OFDM y lo convoluciona con el
% canal
% -----------------------------------------------------------
% function s=TxOFDM(A,N,M,dib)
% -----------------------------------------------------------
% s: símbolo OFDM a transmitir (longitud: N+M)
% h: respuesta impulsiva variante del canal
% i: índice del símbolo OFDM en la señal completa
% -----------------------------------------------------------

function [sf,zf]=CanalLTV(s,h,i)

hi=h(i,:); % Resp imp. i-ésima
sincro=energyTrigger(hi,20);
hi_sincro=hi(sincro:end);
sf=filter(hi_sincro,1,s);