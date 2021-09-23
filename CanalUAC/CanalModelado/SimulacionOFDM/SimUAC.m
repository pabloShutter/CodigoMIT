% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% SimUAC.m
% Script que simula un sistema de comunicaciones digitales
% en el que se transmite una secuencia de s�m5bolos OFDM 
% utilizando codificaci�n QPSK y modulaci�n DMT
% -----------------------------------------------------------
paramSimulacionUAC
filt=1;
noise=0;
N0_2=5e-13;

%% ************ Respuesta del canal UAC ********************
[H,h]=RespuestaUAC(fmin,fmax,d0,wTx,wRx,w,pmax, 1);
ret=length(h)-1;

%% **************  MODULADOR DIGITAL *******************
% Generador de bits
nbitxSimb=Nutil*2;           %n� de bits por simb OFDM
nbit=nbitxSimb*nSymb;        %n� total de bits
nSymbQPSK=Nutil*nSymb;       %n� total de simb QPSK

fprintf('n� s�mbolos OFDM =  %d\n',nSymb)
fprintf('n� s�mbolos QPSK =  %d\n',nSymbQPSK)
fprintf('n� bits simulados = %d\n\n',nbit)

tx=zeros(1,(N+M)*nSymb);     % Inicializamos el vector de la se�al OFDM
bb=zeros(1,nbit);            % Vector de bits que se van a transmitir
AA=zeros(1,nSymbQPSK);       % Vector de s�mbolos QPSK que se van a transmitir

iFc=freqIndex(fc,fs,Df);     %�ndice matlab de frec portadora banda �til
iFcH=freqIndex(-fc,fs,Df);   %�ndice matlab frec portadora banda herm�tica

for kSymb=1:nSymb
   BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb)); % Ventana para el vector de bits
   b=round(rand(1,nbitxSimb));  % Secuencia de bits aleatorios equiprobables
   b(1)=1; b(length(b)/2+1)=1; %Para modulaci�n QPSK diferencial (el primero se toma como ref.)
   bb(BitWin)=b; % Se a�ade la ventana de nbitsSymbOFDM bits al vector de bits a transmitir
   
   A=CodSimbQPSK(b,d);                    % Se codifican los bits en s�mbolos QPSK
   QPSKwin=(kSymb-1)*Nutil+(1:Nutil);  % Ventana para el vector de s�mbolos QPSK
   AA(QPSKwin)=A;      % Se a�aden los Nutil s�mbolos QPSK al vector de 
                       % s�mbolos QPSK a transmitir
   s=TxOFDM(A,N,M,Nutil,iFc,iFcH);    % Se modulan los s�mbolos QPSK en un s�mbolo OFDM
   OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); % Ventana para el vector de s�mbolos OFDM
   tx(OFDMwin)=s;                     % escribe un bloque de se�al
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
    Tn=((N+M)*nSymb+ret)/fs; % Duracion de la se�al de ruido
    tn=0:Tm:Tn-Tm;
end
if noise==1
    n=Ruido(tn,N0_2);        % muestras de ruido AWGN
    yy=sf+n;                % Se�al filtrada con ruido AWG
elseif noise==2
    UWAN=RuidoUWAN(tn,N0_2*500,dib); %Factor 500 para N0_2 media en la banda �til sea aprox 10^-3 W/Hz
    yy=sf+UWAN;
else
    yy=sf;
end

%% *************  DEMODULADOR DIGITAL ******************
QQ=zeros(1,nSymbQPSK);  %S�mbolos QPSK recibidos
AAest=QQ;               %S�mbolos QPSK decididos
bbest=zeros(1,nbit); % 2 bits/symb QPSK

sincro=1;
if filt
    sincro=energyTrigger(h,20);      % �ndice retardo inicial respuesta impulsiva
end
rr=yy(sincro:end); %Desprecio el retardo introducido por el filtro

% Estimaci�n del FEQ
FEQ=1;
if filt
    h_FEQ=h(sincro:sincro+N-1);
    H_FEQ=fft(h_FEQ,N);
    FEQ=1./H_FEQ;
end

for kSymb=1:nSymb
    OFDMwin=(kSymb-1)*(N+M)+(1:(N+M)); %Ventana para recibir el s�mbolo OFDM
    r=rr(OFDMwin);
    
    ro=r((1+M):(N+M));  % Desprecia CP
    Q=RxOFDM(ro,FEQ,Nutil, iFc);  % escribe un bloque de se�al
    
    QPSKwin=(kSymb-1)*Nutil+(1:Nutil);
    QQ(QPSKwin)=Q;
    [aest,best]=DecisorQPSK(Q,d);  % Decisor QPSK y bits
    AAest(QPSKwin)=aest;         % Vector de s�mbolos QPSK decididos
    
    BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb));
    bbest(BitWin)=best;          % Vector de bits decididos
end

%% ****************  PRESTACIONES  *********************
% 1. Probabilidad de error de s�mbolo y de bit simuladas
errsim=abs(AAest-AA);     % secuencia de errores de simbolo
errbit=abs(bbest-bb);     % secuencia de errores de bits
 
PM=nnz(errsim)/nSymbQPSK;     % probabilidad de error de simbolo
Pb=nnz(errbit)/nbit;  % probabilidad de error de bit

% 2. Relaci�n Se�al a Error de Decisi�n
MSE=immse(AA,QQ);
SDR_dB=10*log10(MSE);
if noise
    fprintf('DEP de ruido:                           N0_2 = %g W/Hz\n\n',N0_2)
end

% Energ�a relativa de la respuesta al impulso seg�n N y M
energiah=h*h';
    % N
h_N=h(sincro:sincro+N);
energiah_N=h_N*h_N';
e_N=100*energiah_N/energiah;
    % M
h_M=h(sincro:sincro+M);
energiah_M=h_M*h_M';
e_M=100*energiah_M/energiah;

% 3. Otros par�metros
eta=N/(N+M);        % Eficiencia del s�mbolo OFDM
eta_u=Nutil/(N+M);
Tsimb=Tm*(N+M);     % Periodo de s�mbolo OFDM
Rb=nbitxSimb/Tsimb; % R�gimen binario (b/s)
disp('*********************************************')
fprintf('       N=%d ; M=%d ; Nutil=%d\n',N,M,Nutil)
disp('*********************************************')
fprintf('Relaci�n se�al a distorsi�n:            SDR = %4.2f dB\n',SDR_dB)
fprintf('Probabilidad de error de S�mbolo:       PM = %d\n',PM)
fprintf('Probabilidad de error de Bit:           Pb = %d\n',Pb)
fprintf('Eficiencia de s�mbolo OFDM:             n = %4.2f \n',eta)
fprintf('Eficiencia �til de s�mbolo OFDM:        n_u = %4.2f \n',eta_u)
fprintf('Periodo de s�mbolo OFDM:                Tsimb = %4.2f ms\n',Tsimb*1e3)
fprintf('R�gimen binario:                        Rb = %4.2f Kbps\n',Rb/1e3)
fprintf('Energ�a rel de h (N muestras):          e_N = %4.6f',e_N), disp('%')
fprintf('Energ�a rel de h (M muestras):          e_M = %4.6f',e_M), disp('%')


%% **********  REPRESENTACIONES GR�FICAS  **************-
% FigurasUAC
figure('Name','Comprobacion de recepcion')
subplot(1,2,1), plot(A,'o'), hold on, title('Constelaci�n recibida')
plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(QQ,'.'), hold off, axis([-3 3 -3 3]), grid on
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square

subplot(1,2,2), plot(A,'o'), hold on, title('�ltimo s�mbolo OFDM')
plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(Q,'.'), hold off, axis([-3 3 -3 3]), grid on
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square