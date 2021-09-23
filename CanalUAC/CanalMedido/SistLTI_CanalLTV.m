% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% SistLTI_CanalLTV.m
% Script que simula un sistema de comunicaciones digitales 
% diseñado para trabajar con un canal UAC LTI en el que se 
% transmite una secuencia de símbolos OFDM utilizando 
% codificación QPSK y modulación DMT a través de un canal
% UAC medido y variante con el tiempo (LTV). Este script
% intenta probar la tolerancia de un sistema LTI a un canal
% LTV 
% -----------------------------------------------------------
path(path,'Library')

paramSimUAC_LTV
mod='DQPSK';  %modulación empleada: 'QPSK', '2DPSK' o 'DQPSK'

%% **************  MODULADOR DIGITAL *******************
pos=1;  % Índice de la resp al imp de la que se parte
% Generador de bits
if strcmp(mod,'BPSK') || strcmp(mod,'2DPSK')
    nbitxSimb=Nutil*1;     %nº de bits por simb OFDM
elseif strcmp(mod,'QPSK') || strcmp(mod,'DQPSK')
    nbitxSimb=Nutil*2;
end

nbit=nbitxSimb*nSymb;        %nº total de bits
nSymbQPSK=Nutil*nSymb;       %nº total de simb QPSK

tx=zeros(1,(N+M)*nSymb);     % Inicializamos el vector de la señal OFDM
bb=zeros(1,nbit);            % Vector de bits que se van a transmitir
AA=zeros(1,nSymbQPSK);       % Vector de símbolos QPSK que se van a transmitir

iFc=freqIndex(fc,fs,Df);     %índice matlab de frec portadora banda útil
iFcH=freqIndex(-fc,fs,Df);   %índice matlab frec portadora banda hermítica

for kSymb=1:nSymb
   BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb)); % Ventana para el vector de bits
   b=round(rand(1,nbitxSimb));  % Secuencia de bits aleatorios equiprobables de cada símbolo OFDMx
   
   if strcmp(mod,'BPSK')
       b(1)=1;
       A=CodSimbBPSK(b,d);
   elseif strcmp(mod,'2DPSK')
       b(1)=1;
       A=CodSimb2DPSK(b,d);    % Se codifican los bits en símbolos BPSK
   elseif strcmp(mod,'QPSK')
       b(1)=1; b(length(b)/2+1)=1;  % El primer símbolo QPSK lo mandamos de referencia para DQPSK (1,1)
       A=CodSimbQPSK(b,d);    % Se codifican los bits en símbolos QPSK
   elseif strcmp(mod,'DQPSK')
%        b(1)=1; b(length(b)/2+1)=1;  % El primer símbolo QPSK lo mandamos de referencia para DQPSK (1,1)
%        A=CodSimbDQPSK(b,d);    % Se codifican los bits en símbolos QPSK
        dqpskmod = comm.DQPSKModulator('BitInput',true,'SymbolMapping','Binary');
        dqpskdemod = comm.DQPSKDemodulator('BitOutput',true,'SymbolMapping','Binary');
        A = dqpskmod(b.');
   end
   
   bb(BitWin)=b; % Se añade la ventana de nbitsSymbOFDM bits al vector de bits a transmitir
   
   QPSKwin=(kSymb-1)*Nutil+(1:Nutil);  % Ventana para el vector de símbolos QPSK
   AA(QPSKwin)=A;      % Se añaden los Nutil símbolos QPSK al vector de 
                       % símbolos QPSK a transmitir
   s=TxOFDM(A,N,M,Nutil,iFc,iFcH);    % Se modulan los símbolos QPSK en un símbolo OFDM
   
   if ch==1
       sf=CanalLTV(s,h,kSymb+(pos-1));  % Aplica una resp imp a cada símb OFDM tx
   else
       sf=s;
   end
   
   OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); % Ventana para el vector de símbolos OFDM
   tx(OFDMwin)=sf;                     % escribe un bloque de señal
end

%% ********************  CANAL *************************
ps=mean(abs(tx.^2));

% Ruido aditivo
T=(N+M)*nSymb/fs;      % Duracion de la simulacion
tn=0:Tm:T-Tm;

if tipoRuido==1
    n=Ruido(tn,N0_2);        % muestras de ruido AWGN
    rr=tx+n;                % Señal filtrada con ruido AWG
    pn=mean(abs(n.^2));
elseif tipoRuido==2
    n=RuidoUWAN(tn,N0_2*2,0);
    rr=tx+n;
    pn=mean(abs(n.^2));
else
    rr=tx;
end

%% *************  DEMODULADOR DIGITAL ******************
QQ=zeros(1,nSymbQPSK);  %Símbolos QPSK recibidos
AAest=QQ;               %Símbolos QPSK decididos
bbest=zeros(1,nbit); % 2 bits/symb QPSK

% Estimación del FEQ
hi=h(pos,:); % Resp imp. i-ésima
sincro=energyTrigger(hi,20);
hi_sincro=[ hi(sincro:end) zeros(1,N-length(hi(sincro:end)))];

FEQ=1;
h_FEQ=hi_sincro(1:N);
if ch            
    H_FEQ=fft(h_FEQ,N);
    FEQ=1./H_FEQ; %Igualador ZF
%             FEQ = (conj(H_FEQ)./(abs(H_FEQ).^2 +pn/ps )); %Igualador MMSE
end

errsim=zeros(nSymb,Nutil);
errbit=zeros(nSymb,nbitxSimb);
PM=zeros(1,nSymb);
Pb=zeros(1,nSymb);
SDR_vector=zeros(nSymb,Nutil);
SDR=zeros(1,nSymb);
SDRdB=zeros(1,nSymb);

for kSymb=1:nSymb
    OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el símbolo OFDM
    r=rr(OFDMwin);

    ro=r((1+M):(N+M));  % Desprecia CP
    Q=RxOFDM(ro,FEQ,Nutil, iFc);  % escribe un bloque de señal        

    QPSKwin=(kSymb-1)*Nutil+(1:Nutil);
    QQ(QPSKwin)=Q;
    
    if strcmp(mod,'BPSK')
        [aest,best]=DecisorBPSK(Q,d);        
    elseif strcmp(mod,'2DPSK')
        [aest,best]=Decisor2DPSK(Q,d);  % Decisor QPSK y bits
    elseif strcmp(mod,'QPSK')
        [aest,best]=DecisorQPSK(Q,d);  % Decisor QPSK y bits
    elseif strcmp(mod,'DQPSK')
        best = dqpskdemod(Q.').';
        aest=zeros(1,Nutil);
    end

    AAest(QPSKwin)=aest;         % Vector de símbolos QPSK decididos

    BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb));
    bbest(BitWin)=best;          % Vector de bits decididos
    
    % Calculo y almaceno la PM y la Pb del símbolo k-ésimo
    errsim(kSymb,:)=abs(aest-AA((kSymb-1)*Nutil+(1:Nutil)));
    errbit(kSymb,:)=abs(best-bb((kSymb-1)*nbitxSimb+(1:(nbitxSimb))));
    
    PM(kSymb)=nnz(errsim(kSymb,:))/Nutil;     % probabilidad de error de simbolo
    Pb(kSymb)=nnz(errbit(kSymb,:))/nbitxSimb;          % probabilidad de error de bit
    
%     fprintf('Probabilidad de error de Bit del símbolo #%d :           Pb    = %d\n',kSymb,Pb(kSymb))

    % Calculo y almaceno la SDR del símbolo k-ésimo
    for s=1:Nutil
       SDR_vector(kSymb,s)=norm(AA((kSymb-1)*Nutil+s))/(norm(AA((kSymb-1)*Nutil+s)-Q(s),2));
    end
    SDR(kSymb)=mean(SDR_vector(kSymb,:));
    SDRdB(kSymb)=10*log10(SDR(kSymb));
    
    if dib
       % Represento la constelación del símbolo k-ésimo recibido
        figure('Name',['Símbolo num. ' num2str(kSymb)]), plot(Q,'.'), hold on, title('Último símbolo OFDM')
        plot(A,'o','Linewidth',1.2), hold off, grid on%, axis([-400 400 -400 400])
        xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square 
        pause
    end
end

%% ****************  PRESTACIONES  *********************
% 1. Probabilidad de error de símbolo y de bit simuladas
Pb_mean=mean(Pb);
PM_mean=mean(PM);

% 2. Relación Señal a Error de Decisión por cada símbolo OFDM recibido
SDR_mean=mean(SDR);
SDRdB_mean=mean(SDRdB);

% SNR a la entrada del receptor
P_tx=mean(abs(tx.^2));
P_rx=mean(abs(sf.^2));

% Energía relativa de la respuesta al impulso según N y M
energiah=h(1,:)*h(1,:)';
    % N
energiah_N=h_FEQ*h_FEQ';
e_N=100*energiah_N/energiah;
    % M
h_M=h(1,1:M);
energiah_M=h_M*h_M';
e_M=100*energiah_M/energiah;

e_h=zeros(1,N); %Energía acumulada
for k=1:N
    hh=h(1,1:k);
    energiahh=hh*hh';
    e_h(k)=100*energiahh/energiah;
end

% 3. Otros parámetros
eta=N/(N+M);        % Eficiencia del símbolo OFDM
eta_u=Nutil/(N+M);
Tsimb=Tm*(N+M);     % Periodo de símbolo OFDM
Rb=nbitxSimb/Tsimb; % Régimen binario (b/s)
disp('*********************************************')
fprintf('       N=%d ; M=%d ; Nutil=%d\n',N,M,Nutil)
disp('*********************************************')
%         fprintf('Eficiencia de símbolo OFDM:             n     = %4.2f \n',eta)
%         fprintf('Eficiencia útil de símbolo OFDM:        n_u   = %4.2f \n',eta_u)
        fprintf('Periodo de símbolo OFDM:                Tsimb = %4.2f ms\n',Tsimb*1e3)
%         fprintf('Régimen binario:                        Rb    = %4.2f Kbps\n',Rb/1e3)
%         fprintf('Energía rel de h (N muestras):          e_N   = %4.6f',e_N), disp('%')
%         fprintf('Energía rel de h (M muestras):          e_M   = %4.6f',e_M), disp('%')
        fprintf('Tiempo de coherencia del canal:         Tcoh  = 33.52 ms\n\n')

%% **********  REPRESENTACIONES GRÁFICAS  **************
% close all
% Probabilidad de error de bit por símbolo OFDM recibido
figure('Name','BER por símbolo OFDM recibido')
plot(Pb), grid on, hold on, plot(Pb,'r.'), ylim([0 1])
xlabel('Índice de símbolo OFDM recibido'),title('BER/simbOFDM')

% SDNR por portadora
% figure('Name','SDR por portadora en el primer símbolo OFDM recibido')
% plot(20*log10(SDR_vector(1,:))), grid on
% xlabel('Índice de portadora'),title('SDR/portadora (dB)')

% figure('Name','SDR por símbolo OFDM recibido')
% plot(20*log10(SDR)), grid on
% xlabel('Índice de símbolo OFDM recibido'), title('SDR/símbOFDM (dB)')

%Comprobación de la recepción
figure('Name','Comprobacion de recepcion')
% subplot(1,2,1), 
plot(real(QQ),imag(QQ),'.'), hold on, title('Constelación recibida')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%,l=real(AA(1))*2.5; axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


% figure('Name', 'Último símbolo OFDM recibido'), plot(real(QQ),imag(QQ),'.'), hold on, title('Último símbolo OFDM')
% % plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
% plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%, axis([-l l -l l])
% xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
% ylim([-1 1])