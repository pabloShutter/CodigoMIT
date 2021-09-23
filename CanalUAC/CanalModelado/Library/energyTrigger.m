% -----------------------------------------------------------
% energyTrigger.m
%
% Devuelve el índice a partir del cual comienza la respuesta
% al impulso del canal eliminando el retardo inicial. Se determina
% el índice mediante el cálculo de la energía compredida en una
% ventana móvil de tamaño wSize.
% -----------------------------------------------------------
% function index=delayTrigger(h)
% -----------------------------------------------------------
% Entrada/s:
%   h: Respuesta al impulso completa del canal
%   wSize: Tamaño de la ventana de cálculo de energía
% -----------------------------------------------------------
function index=energyTrigger(h,wSize)

tam=length(h);  %long resp al impulso
e_h=h*h';     %energia de la resp al impulso completa

e=zeros(1,tam-wSize);
for i=1:tam-wSize
    vent=h(i:i+wSize);  %ventana
    e_vent=vent*vent';  %energía de la ventana
    e(i)=e_vent/e_h;    %porcentaje frente a la energía total
end

index=find(e>0.001,1); %Condición: se dispara cuando la energía supera el umbral del 0.5%

% figure, subplot(2,1,1), plot(abs(h)), grid on, title('Respuesta al impulso')
% subplot(2,1,2), plot(e), grid on