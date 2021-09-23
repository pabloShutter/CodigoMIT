% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% SistLTI_CanalLTV.m
% Script que simula un sistema de comunicaciones digitales 
% dise�ado para trabajar con un canal UAC LTI en el que se 
% transmite una secuencia de s�mbolos OFDM utilizando 
% codificaci�n QPSK y modulaci�n DMT a trav�s de un canal
% UAC medido y variante con el tiempo (LTV). Este script
% intenta probar la tolerancia de un sistema LTI a un canal
% LTV 
% -----------------------------------------------------------
path(path,'Library')

paramSimUAC_LTV
mod='DQPSK';  %modulaci�n empleada: 'QPSK', '2DPSK' o 'DQPSK'

%% **************  MODULADOR DIGITAL *******************
pos=1;  % �ndice de la resp al imp de la que se parte
% Generador de bits
if strcmp(mod,'BPSK') || strcmp(mod,'2DPSK')
    nbitxSimb=Nutil*1;     %n� de bits por simb OFDM
elseif strcmp(mod,'QPSK') || strcmp(mod,'DQPSK')
    nbitxSimb=Nutil*2;
end

nbit=nbitxSimb*nSymb;        %n� total de bits
nSymbQPSK=Nutil*nSymb;       %n� total de simb QPSK

tx=zeros(1,(N+M)*nSymb);     % Inicializamos el vector de la se�al OFDM
bb=zeros(1,nbit);            % Vector de bits que se van a transmitir
AA=zeros(1,nSymbQPSK);       % Vector de s�mbolos QPSK que se van a transmitir

iFc=freqIndex(fc,fs,Df);     %�ndice matlab de frec portadora banda �til
iFcH=freqIndex(-fc,fs,Df);   %�ndice matlab frec portadora banda herm�tica

for kSymb=1:nSymb
   BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb)); % Ventana para el vector de bits
   b=round(rand(1,nbitxSimb));  % Secuencia de bits aleatorios equiprobables de cada s�mbolo OFDMx
   
   if strcmp(mod,'BPSK')
       b(1)=1;
       A=CodSimbBPSK(b,d);
   elseif strcmp(mod,'2DPSK')
       b(1)=1;
       A=CodSimb2DPSK(b,d);    % Se codifican los bits en s�mbolos BPSK
   elseif strcmp(mod,'QPSK')
       b(1)=1; b(length(b)/2+1)=1;  % El primer s�mbolo QPSK lo mandamos de referencia para DQPSK (1,1)
       A=CodSimbQPSK(b,d);    % Se codifican los bits en s�mbolos QPSK
   elseif strcmp(mod,'DQPSK')
%        b(1)=1; b(length(b)/2+1)=1;  % El primer s�mbolo QPSK lo mandamos de referencia para DQPSK (1,1)
%        A=CodSimbDQPSK(b,d);    % Se codifican los bits en s�mbolos QPSK
        dqpskmod = comm.DQPSKModulator('BitInput',true,'SymbolMapping','Binary');
        dqpskdemod = comm.DQPSKDemodulator('BitOutput',true,'SymbolMapping','Binary');
        A = dqpskmod(b.');
   end
   
   bb(BitWin)=b; % Se a�ade la ventana de nbitsSymbOFDM bits al vector de bits a transmitir
   
   QPSKwin=(kSymb-1)*Nutil+(1:Nutil);  % Ventana para el vector de s�mbolos QPSK
   AA(QPSKwin)=A;      % Se a�aden los Nutil s�mbolos QPSK al vector de 
                       % s�mbolos QPSK a transmitir
   s=TxOFDM(A,N,M,Nutil,iFc,iFcH);    % Se modulan los s�mbolos QPSK en un s�mbolo OFDM
   
   if ch==1
       sf=CanalLTV(s,h,kSymb+(pos-1));  % Aplica una resp imp a cada s�mb OFDM tx
   else
       sf=s;
   end
   
   OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); % Ventana para el vector de s�mbolos OFDM
   tx(OFDMwin)=sf;                     % escribe un bloque de se�al
end

%% ********************  CANAL *************************
ps=mean(abs(tx.^2));

% Ruido aditivo
T=(N+M)*nSymb/fs;      % Duracion de la simulacion
tn=0:Tm:T-Tm;

if tipoRuido==1
    n=Ruido(tn,N0_2);        % muestras de ruido AWGN
    rr=tx+n;                % Se�al filtrada con ruido AWG
    pn=mean(abs(n.^2));
elseif tipoRuido==2
    n=RuidoUWAN(tn,N0_2*2,0);
    rr=tx+n;
    pn=mean(abs(n.^2));
else
    rr=tx;
end

%% *************  DEMODULADOR DIGITAL ******************
QQ=zeros(1,nSymbQPSK);  %S�mbolos QPSK recibidos
AAest=QQ;               %S�mbolos QPSK decididos
bbest=zeros(1,nbit); % 2 bits/symb QPSK

% Estimaci�n del FEQ
hi=h(pos,:); % Resp imp. i-�sima
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
    OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el s�mbolo OFDM
    r=rr(OFDMwin);

    ro=r((1+M):(N+M));  % Desprecia CP
    Q=RxOFDM(ro,FEQ,Nutil, iFc);  % escribe un bloque de se�al        

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

    AAest(QPSKwin)=aest;         % Vector de s�mbolos QPSK decididos

    BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb));
    bbest(BitWin)=best;          % Vector de bits decididos
    
    % Calculo y almaceno la PM y la Pb del s�mbolo k-�simo
    errsim(kSymb,:)=abs(aest-AA((kSymb-1)*Nutil+(1:Nutil)));
    errbit(kSymb,:)=abs(best-bb((kSymb-1)*nbitxSimb+(1:(nbitxSimb))));
    
    PM(kSymb)=nnz(errsim(kSymb,:))/Nutil;     % probabilidad de error de simbolo
    Pb(kSymb)=nnz(errbit(kSymb,:))/nbitxSimb;          % probabilidad de error de bit
    
%     fprintf('Probabilidad de error de Bit del s�mbolo #%d :           Pb    = %d\n',kSymb,Pb(kSymb))

    % Calculo y almaceno la SDR del s�mbolo k-�simo
    for s=1:Nutil
       SDR_vector(kSymb,s)=norm(AA((kSymb-1)*Nutil+s))/(norm(AA((kSymb-1)*Nutil+s)-Q(s),2));
    end
    SDR(kSymb)=mean(SDR_vector(kSymb,:));
    SDRdB(kSymb)=10*log10(SDR(kSymb));
    
    if dib
       % Represento la constelaci�n del s�mbolo k-�simo recibido
        figure('Name',['S�mbolo num. ' num2str(kSymb)]), plot(Q,'.'), hold on, title('�ltimo s�mbolo OFDM')
        plot(A,'o','Linewidth',1.2), hold off, grid on%, axis([-400 400 -400 400])
        xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square 
        pause
    end
end

%% ****************  PRESTACIONES  *********************
% 1. Probabilidad de error de s�mbolo y de bit simuladas
Pb_mean=mean(Pb);
PM_mean=mean(PM);

% 2. Relaci�n Se�al a Error de Decisi�n por cada s�mbolo OFDM recibido
SDR_mean=mean(SDR);
SDRdB_mean=mean(SDRdB);

% SNR a la entrada del receptor
P_tx=mean(abs(tx.^2));
P_rx=mean(abs(sf.^2));

% Energ�a relativa de la respuesta al impulso seg�n N y M
energiah=h(1,:)*h(1,:)';
    % N
energiah_N=h_FEQ*h_FEQ';
e_N=100*energiah_N/energiah;
    % M
h_M=h(1,1:M);
energiah_M=h_M*h_M';
e_M=100*energiah_M/energiah;

e_h=zeros(1,N); %Energ�a acumulada
for k=1:N
    hh=h(1,1:k);
    energiahh=hh*hh';
    e_h(k)=100*energiahh/energiah;
end

% 3. Otros par�metros
eta=N/(N+M);        % Eficiencia del s�mbolo OFDM
eta_u=Nutil/(N+M);
Tsimb=Tm*(N+M);     % Periodo de s�mbolo OFDM
Rb=nbitxSimb/Tsimb; % R�gimen binario (b/s)
disp('*********************************************')
fprintf('       N=%d ; M=%d ; Nutil=%d\n',N,M,Nutil)
disp('*********************************************')
%         fprintf('Eficiencia de s�mbolo OFDM:             n     = %4.2f \n',eta)
%         fprintf('Eficiencia �til de s�mbolo OFDM:        n_u   = %4.2f \n',eta_u)
        fprintf('Periodo de s�mbolo OFDM:                Tsimb = %4.2f ms\n',Tsimb*1e3)
%         fprintf('R�gimen binario:                        Rb    = %4.2f Kbps\n',Rb/1e3)
%         fprintf('Energ�a rel de h (N muestras):          e_N   = %4.6f',e_N), disp('%')
%         fprintf('Energ�a rel de h (M muestras):          e_M   = %4.6f',e_M), disp('%')
        fprintf('Tiempo de coherencia del canal:         Tcoh  = 33.52 ms\n\n')

%% **********  REPRESENTACIONES GR�FICAS  **************
% close all
% Probabilidad de error de bit por s�mbolo OFDM recibido
figure('Name','BER por s�mbolo OFDM recibido')
plot(Pb), grid on, hold on, plot(Pb,'r.'), ylim([0 1])
xlabel('�ndice de s�mbolo OFDM recibido'),title('BER/simbOFDM')

% SDNR por portadora
% figure('Name','SDR por portadora en el primer s�mbolo OFDM recibido')
% plot(20*log10(SDR_vector(1,:))), grid on
% xlabel('�ndice de portadora'),title('SDR/portadora (dB)')

% figure('Name','SDR por s�mbolo OFDM recibido')
% plot(20*log10(SDR)), grid on
% xlabel('�ndice de s�mbolo OFDM recibido'), title('SDR/s�mbOFDM (dB)')

%Comprobaci�n de la recepci�n
figure('Name','Comprobacion de recepcion')
% subplot(1,2,1), 
plot(real(QQ),imag(QQ),'.'), hold on, title('Constelaci�n recibida')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%,l=real(AA(1))*2.5; axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


% figure('Name', '�ltimo s�mbolo OFDM recibido'), plot(real(QQ),imag(QQ),'.'), hold on, title('�ltimo s�mbolo OFDM')
% % plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
% plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%, axis([-l l -l l])
% xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
% ylim([-1 1])