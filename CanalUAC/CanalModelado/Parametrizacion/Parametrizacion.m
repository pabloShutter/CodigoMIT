% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% Parametrizacion.m
% Script que simula un sistema de comunicaciones digitales
% en el que se transmite una secuencia de símbolos OFDM 
% utilizando codificación QPSK y modulación DMT con la idea
% de conseguir un diseño óptimo
% -----------------------------------------------------------

paramSimulacionUAC

%% ************ Respuesta del canal UAC ********************
[H,h]=RespuestaUAC(fmin,fmax,d0,wTx,wRx,w,pmax, dib);
ret=length(h)-1;

%% **************  MODULADOR DIGITAL *******************
% Generador de bits
nbitxSimb=Nutil*2;           %nº de bits por simb OFDM
nbit=nbitxSimb*nSymb;        %nº total de bits
nSymbQPSK=Nutil*nSymb;       %nº total de simb QPSK

fprintf('nº símbolos OFDM =  %d\n',nSymb)
fprintf('nº símbolos QPSK =  %d\n',nSymbQPSK)
fprintf('nº bits simulados = %d\n\n',nbit)

tx=zeros(1,(N+M)*nSymb);     % Inicializamos el vector de la señal OFDM
bb=zeros(1,nbit);            % Vector de bits que se van a transmitir
AA=zeros(1,nSymbQPSK);       % Vector de símbolos QPSK que se van a transmitir

iFc=freqIndex(fc,fs,Df);     %índice matlab de frec portadora banda útil
iFcH=freqIndex(-fc,fs,Df);   %índice matlab frec portadora banda hermítica

for kSymb=1:nSymb
   BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb)); % Ventana para el vector de bits
   b=round(rand(1,nbitxSimb));  % Secuencia de bits aleatorios equiprobables
   bb(BitWin)=b; % Se añade la ventana de nbitsSymbOFDM bits al vector de bits a transmitir
   
   A=CodSimbQPSK(b,d);                    % Se codifican los bits en símbolos QPSK
   QPSKwin=(kSymb-1)*Nutil+(1:Nutil);  % Ventana para el vector de símbolos QPSK
   AA(QPSKwin)=A;      % Se añaden los Nutil símbolos QPSK al vector de 
                       % símbolos QPSK a transmitir
   s=TxOFDM(A,N,M,Nutil,iFc,iFcH);    % Se modulan los símbolos QPSK en un símbolo OFDM
   OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); % Ventana para el vector de símbolos OFDM
   tx(OFDMwin)=s;                     % escribe un bloque de señal
end

%% ********************  CANAL *************************
% Filtrado respuesta al impulso h[n]
if filt
    sf=conv(h,tx);
else
    sf=tx;
end

% Ruido aditivo
tn=t;
if filt
    Tn=((N+M)*nSymb+ret)/fs; % Duracion de la señal de ruido
    tn=0:Tm:Tn-Tm;
end
if noise==1
    n=Ruido(tn,N0_2);        % muestras de ruido AWGN
    yy=sf+n;                % Señal filtrada con ruido AWG
    pn=mean(abs(n.^2));
elseif noise==2
    n=RuidoUWAN(tn,N0_2*2,dib);
    yy=sf+n;
    pn=mean(abs(n.^2));
else
    yy=sf;
end
%% *************  DEMODULADOR DIGITAL ******************
QQ=zeros(1,nSymbQPSK);  %Símbolos QPSK recibidos
AAest=QQ;               %Símbolos QPSK decididos
bbest=zeros(1,nbit); % 2 bits/symb QPSK

sincro=1;
if filt
    sincro=energyTrigger(h,20);      % índice retardo inicial respuesta impulsiva
end
rr=yy(sincro:end); %Desprecio el retardo introducido por el filtro

% Estimación del FEQ
FEQ=1;
if filt
    h_FEQ=h(sincro:end);
    H_FEQ=fft(h_FEQ,N);
    FEQ=1./H_FEQ;
%     FEQ = (H_FEQ./(abs(H_FEQ).^2 + (pn/ps)))';
%     FEQ = (conj(H_FEQ)./(abs(H_FEQ).^2 + pn/ps )); %Igualador MMSE
end

for kSymb=1:nSymb
    OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el símbolo OFDM
    r=rr(OFDMwin);
    
    ro=r((1+M):(N+M));  % Desprecia CP
    Q=RxOFDM(ro,FEQ,Nutil, iFc);  % escribe un bloque de señal
    
    QPSKwin=(kSymb-1)*Nutil+(1:Nutil);
    QQ(QPSKwin)=Q;
    [aest,best]=DecisorQPSK(Q,d);  % Decisor QPSK y bits
    AAest(QPSKwin)=aest;         % Vector de símbolos QPSK decididos
    
    BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb));
    bbest(BitWin)=best;          % Vector de bits decididos
end

%% ****************  PRESTACIONES  *********************
% 1. Probabilidad de error de símbolo y de bit simuladas
errsim=abs(AAest-AA);     % secuencia de errores de simbolo
errbit=abs(bbest-bb);     % secuencia de errores de bits
 
PM=nnz(errsim)/nSymbQPSK;     % probabilidad de error de simbolo
Pb=nnz(errbit)/nbit;  % probabilidad de error de bit

% 2. Relación Señal a Error de Decisión
% MSE2=immse(AA,QQ);
SDR_vector=zeros(1,nSymbQPSK);
for s=1:nSymbQPSK
   SDR_vector(s)=norm(AA(s))/(norm(AA(s)-QQ(s),2));
end
SDR=mean(SDR_vector);
SDR_dB=10*log10(SDR);

if noise
    fprintf('DEP de ruido:                           N0_2 = %g W/Hz\n\n',N0_2)
end

% SNR a la entrada del receptor
P_tx=mean(abs(tx.^2));
P_rx=mean(abs(sf.^2));
if noise
    Pn_rx=mean(abs(n.^2));
    SNR_rx=(P_rx*N)/(Pn_rx*Nutil);
    SNR_rx_dB=10*log10(SNR_rx);
end

% Energía relativa de la respuesta al impulso según N y M
energiah=h*h';
    % N
h_N=h(sincro:sincro+N);
energiah_N=h_N*h_N';
e_N=100*energiah_N/energiah;
    % M
h_M=h(sincro:sincro+M);
energiah_M=h_M*h_M';
e_M=100*energiah_M/energiah;

e_h=zeros(1,N); %Energía acumulada
for k=1:N
    hh=h(sincro:sincro+k);
    energiahh=hh*hh';
    e_h(k)=100*energiahh/energiah;
end

% figure('Name','Energía de h en función de M')
% plot(e_h), grid on, xlabel('Long. CP (M)'), ylabel('energía de h (%)')

% 3. Otros parámetros
eta=N/(N+M);        % Eficiencia del símbolo OFDM
eta_u=Nutil/(N+M);
Tsimb=Tm*(N+M);     % Periodo de símbolo OFDM
Rb=nbitxSimb/Tsimb; % Régimen binario (b/s)
disp('*********************************************')
fprintf('       N=%d ; M=%d ; Nutil=%d\n',N,M,Nutil)
disp('*********************************************')
fprintf('Relación señal a distorsión:            SDR   = %4.2f dB\n',SDR_dB)
if noise 
    fprintf('Relación señal a ruido entrada Rx:      SNRrx = %4.2f dB\n',SNR_rx_dB)
end
fprintf('Probabilidad de error de Símbolo:       PM    = %d\n',PM)
fprintf('Probabilidad de error de Bit:           Pb    = %d\n',Pb)
fprintf('Eficiencia de símbolo OFDM:             n     = %4.2f \n',eta)
fprintf('Eficiencia útil de símbolo OFDM:        n_u   = %4.2f \n',eta_u)
fprintf('Periodo de símbolo OFDM:                Tsimb = %4.2f ms\n',Tsimb*1e3)
fprintf('Régimen binario:                        Rb    = %4.2f Kbps\n',Rb/1e3)
fprintf('Energía rel de h (N muestras):          e_N   = %4.6f',e_N), disp('%')
fprintf('Energía rel de h (M muestras):          e_M   = %4.6f',e_M), disp('%')


%% **********  REPRESENTACIONES GRÁFICAS  **************-
% FigurasUAC
% SDNR por portadora 
figure('Name','SDNR por portadora'), plot(1:Nutil, 10*log10(SDR_vector(1:Nutil))), grid on
x=1:Nutil; y=10*log10(SDR_vector(1:Nutil));
p=polyfit(x,y,2);
y1=polyval(p,x);
hold on, plot(x,y1,'r--'), hold off
xlabel('Símbolos QPSK (portadoras)'), ylabel('SNDR (dB)'), title('SNDR/portadora con AWGN')

%Comprobación de la recepción
figure('Name','Comprobacion de recepcion')
subplot(1,2,1), plot(QQ,'.'), hold on, title('Constelación recibida')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(A,'o'), hold off, grid on%,l=real(AA(1))*2.5; axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square

subplot(1,2,2), plot(Q,'.'), hold on, title('Último símbolo OFDM')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(A,'o'), hold off, grid on%, axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square