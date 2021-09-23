% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% SimSecuenciaOFDM.m
% Script que simula un sistema de comunicaciones digitales
% utilizando codificación QPSK y modulación con una secuencia
% de símbolos OFDM
% -----------------------------------------------------------

clear,clc%,close all;
path(path,'Library')
paramSimulacion

nbitxSimb=Nutil*2;           %nº de bits por simb OFDM
nbit=nbitxSimb*nSymb;        %nº total de bits
nSymbQPSK=Nutil*nSymb;       %nº total de simb QPSK

SNRdB=2:2:12;

for i=1:length(SNRdB)
    %% **************  MODULADOR DIGITAL *******************
    % Generador de bits
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
    
    %Debería estar calculado en el Rx pero lo necesito aquí
    FEQ=1;
    if filt
        FEQ=1./(H.');
    end

    SNR=10^(SNRdB(i)/10);
    senyal=FEQyQuitaCP(sf,N,M,Nutil,FEQ,nSymb,iFc,iFcH);
    Er=sum(abs(senyal).^2)*Tm/nSymbQPSK;
    Eb=Er/2;
    N0_2(i)=Eb/(2*SNR);

    % Ruido aditivo
    tn=t;
    if filt
        Tn=((N+M)*nSymb+ret)/fs; % Duracion de la señal de ruido
        tn=0:Tm:Tn-Tm;
    end
    if noise==1
        n=Ruido(tn,N0_2(i));
        rr=sf+n;                % Señal filtrada con ruido AWG
    elseif noise==2
        n=RuidoUWAN(tn,N0_2,1);  % muestras UWAN
        rr=sf+n;                % Señal filtrada con ruido UWA
    else
        rr=sf;
    end

    if noise~=0
        ruido=FEQyQuitaCP(n,N,M,Nutil,1,nSymb,iFc,iFcH);
        ruidoDEP=EstimaDEP(ruido,N,fs,'b');
        nDEP=[ruidoDEP(iFc-Nutil/2:iFc-1) ruidoDEP(iFc+1:iFc+Nutil/2)];
        N0_22=mean(nDEP);
        SNRdB2(i)=10*log10(Eb/(2*N0_22));
    end

    %% *************  DEMODULADOR DIGITAL ******************
    QQ=zeros(1,nSymbQPSK);  %Símbolos QPSK recibidos
    AAest=QQ;               %Símbolos QPSK decididos
    bbest=zeros(1,nbit); % 2 bits/symb QPSK

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
    Pb(i)=nnz(errbit)/nbit;  % probabilidad de error de bit

    % 2. Relación Señal a Error de Decisión
    MSE=immse(AA,QQ);
    SDER_dB=10*log10(MSE);
    
    SNR_vect=zeros(1,nSymbQPSK); %En el canal UAC lo llamo SNDR, pero aquí no hay distorsión, así que se puede considerar solo SNR
    SNR_vect_mean=zeros(1,Nutil);
    for s=1:nSymbQPSK
       SNR_vect(s)=norm(AA(s))/(norm(AA(s)-QQ(s),2));
       SNR_vect_mean(mod(s-1,Nutil)+1)=SNR_vect_mean(mod(s-1,Nutil)+1) + norm(AA(s))/(norm(AA(s)-QQ(s),2));
    end
    SNR_vect_mean=SNR_vect_mean./nSymb;
    SDR=mean(SNR_vect);
    SDR_dB=10*log10(SDR);
end

% Información
if noise
    for i=1:length(SNRdB2)
        fprintf('Relación señal a Ruido:                 SNR = %4.2f dB\n',SNRdB2(i))
        fprintf('Probabilidad de error de Bit:           Pb = %1.3e\n\n',Pb(i))
    end
end
%% **********  REPRESENTACIONES GRÁFICAS  **************-
% 3. Cálculo probabilidad de error de bit teórica
if noise
    EbNoVec=(0:0.1:12.2)';
    berTheory=qfunc(sqrt(2*10.^(EbNoVec/10)));
    figure, semilogy(EbNoVec,berTheory,'b--','linewidth',1), grid on    
    xlabel('Eb/No (dB)'), ylabel('Bit Error Rate (BER)')

    hold on, semilogy(SNRdB2,Pb,'r*', 'linewidth',1)
    legend('Teórica','Canal FIR con FEQ')
end

% Comprobación de la recepción
figure('Name','Comprobacion de recepcion')
plot(QQ,'.'), hold on
plot(AA,'o','Linewidth',2), hold off, grid on
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';

figure, line((1:length(SNR_vect_mean)),10*log10(SNR_vect_mean),'Linewidth',1.5), grid on
axis square, xlabel('Índice de portadora'), ylabel('SNR (dB)')