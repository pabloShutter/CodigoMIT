% -----------------------------------------------------------
% Autor: Pablo Cobacho
% -----------------------------------------------------------
% SIMOFDM_CM.m
% Script que simula un sistema de comunicaciones digitales
% en el que se transmite una secuencia de s�mbolos OFDM 
% utilizando codificaci�n QPSK y modulaci�n DMT a trav�s de
% un canal UAC est�tico medido
% -----------------------------------------------------------

path(path,'Library')
paramSimulacionUAC

%% **************  MODULADOR DIGITAL *******************
N=[4096 8192];
Nutil=[600 1200];
% N=4096; Nutil=600;

for i=1:length(N)
    M=[N(i)*3/8 N(i)/2 N(i)*3/4 N(i)];
%     M=N*3/4;
    for j=1:length(M)
        Df=fs/N(i);
        % Generador de bits
        nbitxSimb=Nutil(i)*2;           %n� de bits por simb OFDM
        nbit=nbitxSimb*nSymb;        %n� total de bits
        nSymbQPSK=Nutil(i)*nSymb;       %n� total de simb QPSK

        tx=zeros(1,(N(i)+M(j))*nSymb);     % Inicializamos el vector de la se�al OFDM
        bb=zeros(1,nbit);            % Vector de bits que se van a transmitir
        AA=zeros(1,nSymbQPSK);       % Vector de s�mbolos QPSK que se van a transmitir

        iFc=freqIndex(fc,fs,Df);     %�ndice matlab de frec portadora banda �til
        iFcH=freqIndex(-fc,fs,Df);   %�ndice matlab frec portadora banda herm�tica

        for kSymb=1:nSymb
           BitWin=(kSymb-1)*nbitxSimb+(1:(nbitxSimb)); % Ventana para el vector de bits
           b=round(rand(1,nbitxSimb));  % Secuencia de bits aleatorios equiprobables
%            b(1)=1; b(length(b)/2+1)=1;  % El primer s�mbolo QPSK lo mandamos de referencia para DQPSK (1,1)
%            b(1)=1;
           bb(BitWin)=b; % Se a�ade la ventana de nbitsSymbOFDM bits al vector de bits a transmitir           
           
           A=CodSimbQPSK(b,d);                    % Se codifican los bits en s�mbolos QPSK
           
           QPSKwin=(kSymb-1)*Nutil(i)+(1:Nutil(i));  % Ventana para el vector de s�mbolos QPSK
           AA(QPSKwin)=A;      % Se a�aden los Nutil s�mbolos QPSK al vector de 
                               % s�mbolos QPSK a transmitir
           s=TxOFDM(A,N(i),M(j),Nutil(i),iFc,iFcH);    % Se modulan los s�mbolos QPSK en un s�mbolo OFDM
           OFDMwin=(kSymb-1)*(N(i)+M(j))+(1:(N(i)+M(j))); % Ventana para el vector de s�mbolos OFDM
           tx(OFDMwin)=s;                     % escribe un bloque de se�al
        end

        %% ********************  CANAL *************************
        % Filtrado respuesta al impulso h[n]
        if ch
            sf=conv(h,tx);
        else
            sf=tx;
        end
        ps=mean(abs(sf.^2));

        % Ruido aditivo
        T=(N(i)+M(j))*nSymb/fs;      % Duracion de la simulacion
        tn=0:Tm:T-Tm;
        if ch
            Tn=((N(i)+M(j))*nSymb+ret)/fs; % Duracion de la se�al de ruido
            tn=0:Tm:Tn-Tm;
        end
        if noise==1
            n=Ruido(tn,N0_2);        % muestras de ruido AWGN
            rr=sf+n;                % Se�al filtrada con ruido AWG
            pn=mean(abs(n.^2));
        elseif noise==2
            n=RuidoUWAN(tn,N0_2*2,0);
            rr=sf+n;
            pn=mean(abs(n.^2));
        else
            rr=sf;
        end
        %% *************  DEMODULADOR DIGITAL ******************
        QQ=zeros(1,nSymbQPSK);  %S�mbolos QPSK recibidos
        AAest=QQ;               %S�mbolos QPSK decididos
        bbest=zeros(1,nbit); % 2 bits/symb QPSK

        sincro=1;

        % Estimaci�n del FEQ
        FEQ=1;
        
        if ch
            h_FEQ=h(1:N(i));
            H_FEQ=fft(h_FEQ,N(i));
            FEQ=1./H_FEQ; %Igualador ZF
%             FEQ = (conj(H_FEQ)./(abs(H_FEQ).^2 +pn/ps )); %Igualador MMSE
        end

        for kSymb=1:nSymb
            OFDMwin=(kSymb-1)*(N(i)+M(j))+(1:(N(i)+M(j))); %Ventana para recibir el s�mbolo OFDM
            r=rr(OFDMwin);

            ro=r((1+M(j)):(N(i)+M(j)));  % Desprecia CP
            Q=RxOFDM(ro,FEQ,Nutil(i), iFc);  % escribe un bloque de se�al

            QPSKwin=(kSymb-1)*Nutil(i)+(1:Nutil(i));
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
        % MSE2=immse(AA,QQ);
        SDR_vector=zeros(1,nSymbQPSK);
        for s=1:nSymbQPSK
           SDR_vector(s)=norm(AA(s))/(norm(AA(s)-QQ(s),2));
        end
        SDR=mean(SDR_vector);
        SDR_dB=10*log10(SDR);

%         if noise
%             fprintf('DEP de ruido:                           N0_2 = %g W/Hz\n\n',N0_2)
%         end

        % SNR a la entrada del receptor
        P_tx=mean(abs(tx.^2));
        P_rx=mean(abs(sf(sincro:sincro+(nSymb*(N(i)+M(j)))).^2));
        if noise
            Pn_rx=mean(abs(n.^2));
            SNRp_rx=(P_rx*N(i))/(Pn_rx*Nutil(i)); %SNR/portadora a la entrada del rx
            SNRp_rx_dB=10*log10(SNRp_rx);
            SNR_rx=(P_rx)/(Pn_rx);           %SNR a la entrada del rx
            SNR_rx_dB=10*log10(SNR_rx);
        end

        % Energ�a relativa de la respuesta al impulso seg�n N y M
        if ch==1
            energiah=h*h';
                % N
            energiah_N=h_FEQ*h_FEQ';
            e_N=100*energiah_N/energiah;
                % M
            h_M=h(1:M(j));
            energiah_M=h_M*h_M';
            e_M=100*energiah_M/energiah;

            e_h=zeros(1,N(i)); %Energ�a acumulada
            for k=1:N(i)
                hh=h(1:k);
                energiahh=hh*hh';
                e_h(k)=100*energiahh/energiah;
            end
        end
        
        % 3. Otros par�metros
        eta=N(i)/(N(i)+M(j));        % Eficiencia del s�mbolo OFDM
        eta_u=Nutil(i)/(N(i)+M(j));
        Tsimb=Tm*(N(i)+M(j));     % Periodo de s�mbolo OFDM
        Rb=nbitxSimb/Tsimb; % R�gimen binario (b/s)
        disp('*********************************************')
        fprintf('       N=%d ; M=%d ; Nutil=%d\n',N(i),M(j),Nutil(i))
        disp('*********************************************')
        fprintf('Relaci�n se�al a distorsi�n:            SDR   = %4.1f dB\n',SDR_dB)
        if noise 
            fprintf('SNR/portadora a la entrada del Rx:      SNRp_rx = %4.1f dB\n',SNRp_rx_dB)
            fprintf('SNR total a la entrada del Rx:          SNR_rx = %4.1f dB\n',SNR_rx_dB)
        end
%         fprintf('Probabilidad de error de S�mbolo:       PM    = %d\n',PM)
        fprintf('Probabilidad de error de Bit:           Pb    = %.2e\n',Pb)
%         fprintf('Eficiencia de s�mbolo OFDM:             n     = %4.2f \n',eta)
%         fprintf('Eficiencia �til de s�mbolo OFDM:        n_u   = %4.2f \n',eta_u)
%         fprintf('Periodo de s�mbolo OFDM:                Tsimb = %4.2f ms\n',Tsimb*1e3)
%         fprintf('R�gimen binario:                        Rb    = %4.2f Kbps\n',Rb/1e3)
%         if ch==1
%             fprintf('Energ�a rel de h (N muestras):          e_N   = %4.6f',e_N), disp('%')
%             fprintf('Energ�a rel de h (M muestras):          e_M   = %4.6f',e_M), disp('%')
%         end
    end
end

%% **********  REPRESENTACIONES GR�FICAS  **************-
%Energ�a de la respuesta al impulso en funci�n de M
% figure('Name','Energ�a de h en funci�n de M')
% plot(e_h), grid on, xlabel('Long. CP (M)'), ylabel('energ�a de h (%)')

% FigurasUAC
% SDNR por portadora 
% figure('Name','SDNR por portadora'), plot(1:Nutil, 10*log10(SDR_vector(1:Nutil))), grid on
% x=1:Nutil; y=10*log10(SDR_vector(1:Nutil));
% p=polyfit(x,y,2);
% y1=polyval(p,x);
% hold on, plot(x,y1,'r--'), hold off
% xlabel('S�mbolos QPSK (portadoras)'), ylabel('SNDR (dB)'), title('SNDR/portadora con AWGN')

%Comprobaci�n de la recepci�n
figure('Name','Comprobacion de recepcion')
subplot(1,2,1), plot(real(QQ),imag(QQ),'.'), hold on, title('Constelaci�n recibida')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%,l=real(AA(1))*2.5; axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% ylim([-0.5 0.5])

subplot(1,2,2), plot(real(QQ),imag(QQ),'.'), hold on, title('�ltimo s�mbolo OFDM')
% plot(-3:.1:3, zeros(1,61),'r--', zeros(1,61),-3:.1:3,'r--','LineWidth',2)
plot(real(AA),imag(AA),'o','Linewidth',1.2), hold off, grid on%, axis([-l l -l l])
xlabel('Real\{a_k\}'), ylabel('Imag\{a_k\}'), axis square
% ylim([-1 1])