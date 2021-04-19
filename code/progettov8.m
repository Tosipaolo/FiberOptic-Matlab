%% Parametri

ptx_bit1_dbm = 0; %potenza trasmessa sul bit 1 in dBm
ptx_bit1 = 10^(-3)*10^(ptx_bit1_dbm/10);   %conversione in watt
epsilon = 9.7; %rapporto di estinzione lineare approssimato dal valore 9.69
ptx_bit0 = ptx_bit1 / epsilon; 
Lcanale_km = 74;  %lughezza canale definita in km 
LEN= 2048;     %numero di bit
dt = 1e-12;    %risoluzione temporale
tsimbolo= 100e-12; %tempo di bit
tempo_totale = LEN*tsimbolo; %ampiezza dell'asse temporale
numerocampioni = (tempo_totale/dt)+1; %numero di campioni 


Bn= 7.5e9; %banda equivalente del rumore
responsivity= 1; % indice conversione [A/W]
NEC= 20e-12; % [A/sqrt(Hz)] Noise Equivality current
f_portottica=2e14;
N=tempo_totale/tsimbolo; %numero bit

frequenza_totale = 1/dt; %ampiezza asse delle frequenze [Hz]
df= (1/tempo_totale); %passo in frequenza o risoluzione frequenziale [Hz]

f_s=[ 0  : df : frequenza_totale];
f_simmetrico = [-(frequenza_totale / 2) : df : frequenza_totale / 2 ];

bit_time = tempo_totale/LEN;%secondi
bit_intervalli = bit_time/dt; %intervalli per bit

t_s = 0:dt:tempo_totale; %array tempo-continuo

%% sequenza di bit casuali 
bit = single.empty(LEN,0);

bit = randi([0 1], 1,LEN);

dig_sig = zeros(size(t_s));

%% creazione segnale digitale dai valori della serie di bit

for i = 1:LEN
    if(bit(i)) == 1
         dig_sig(((i-1)*bit_intervalli) +1 : i*bit_intervalli) = bit(i);
    else
         dig_sig(((i-1)*bit_intervalli) +1 : i*bit_intervalli) = 1/3;  
    end

end

if bit(LEN) == 1
    dig_sig(tempo_totale/dt + 1) = bit(LEN);
else
    dig_sig(tempo_totale/dt +1) = 1/3;
end

figure (1);   %plotta il segnale rettangolare ideale
plot(t_s, dig_sig);
title('Segnale Rettangolare','color','blue');
grid on;

%% creazione coseno rialzato

rcos = zeros(size(t_s));
singlercos = zeros(size(t_s));
alpha=1;
rcos = ((1/2)* (1+cos(2*pi*t_s/(alpha*tsimbolo))))/50  ;


for i=51:151
    singlercos(i) = rcos(i);
end

% figure (2);  %singolo coseno rialzato da 0.05 a 0.15
% plot(t_s, singlercos);
% grid on; 


%% convoluzione discreta nel dominio del tempo/ segnale trasmesso

DIG_SIG= fft(dig_sig);        %fft del mio segnale rettang.
RCOS= fft(singlercos);        %fft del mio coseno rialzato singolo
TX= DIG_SIG .* RCOS;       %prodotto nel d.spettrale delle frequenze

tx=ifft(TX);   %fft inversa per riporta nel dominio temporale(segnale tx.)

figure(3);    %segnale trasmesso sul canale (smussato) con traslazione
plot(t_s,tx);
grid on;
title('Segnale smussato non ancora traslato ','color','blue');

%% traslazione per riaccostare il segnale 

tx_shifted = zeros(size(t_s)); 
j=1;
for i= 101 : numerocampioni 
    tx_shifted(j) = tx(i);
    j= j+1;
end

for k= 1 : 100
    tx_shifted(j) = tx(k);
    j=j+1;
end

TX_s = fft(tx_shifted);

figure(4); %plotto il segnale smussato e traslato
plot(t_s, tx_shifted);
grid on;
title('Segnale trasmesso smussatoe tralsato','color','blue');


%% Definzione della risposta in frequenza del CANALE
lungonda = 1550e-9; %lunghezza d'onda in [m]
c= 3e8; %velocità della luce in [m/s]
eta_eff = 1.46; %indice di rifrazione del silicio
D= 16.7e-6; %coefficiente di dispersione cromatica definito in [s/(m^2)]
beta2 = -D*((lungonda^2)*2*pi/c) ; %derivata seconda della risposta in fase
alpha_dbkm = 0.2;  %coefficiente di attenuazione [dB/Km]
alpha = (1/4.343)* alpha_dbkm;  %coefficiente attenuazione in [1/Km]

A_f = exp((-1/2)*alpha*Lcanale_km); %attenuazione che considero a livello di potenza ma non relativamente al canale


B_f = zeros(size(f_simmetrico));  
freq_diappoggio = f_simmetrico ;  %vettore copia f_simmetrico
% B_f(k)= 0.5* beta2.*(((freq_diappoggio(k))^2).*Lcanale_km * 10^3);
for k= 1: length(f_simmetrico);
    B_f(k)= beta2 *((freq_diappoggio(k)).^2)* Lcanale_km * 10^3 ;
end

%% Compensazione con fibra speciale

alpha_specialdbkm = 0.5/6 %coef.di attenuazione in dB/Km
alpha_slineare = alpha_specialdbkm*(1/ 4.343);
alpha_totale = alpha_slineare + alpha;
A_fspecial= exp((-1/2)*alpha_totale*Lcanale_km);

%B_f= zeros(size(f_simmetrico));

%%

Hc_f0 = complex(cos(B_f) , sin(B_f));  %ris.freq.canale(complessa) equivalente in b.base A(f)^2 verrà considerato successivamente al canale come attenuazione di potenza riga 165
%Hc_f1 = ((1+sign(freq_diappoggio+f_portottica))/2).*Hc_f0;  %%inviluppo complesso in b.base elimino fnegative e traslo
Hc_f= fftshift(Hc_f0); %shift della risposta in banda monolatera

figure (13);   %plotto la risposta in frequenza del canale traslata
plot(f_s, Hc_f);
title('risposta in frequenza canale','color','blue');
grid on;
 
RX = TX_s .* Hc_f;  %filtraggio segnale attraverso il canale 
rx= ifft(RX);  %segnale post canale riportato nel dom. temporale

rx_intensity= (abs(rx).^2);  %intesità luminosa 

figure(5);
plot(t_s, rx_intensity); %plotto l'andamento del modulo quadro del segnale post canale
grid on;
title('Segnale dopo il canale ','color','blue');

%% Fase di ricezione, introduzione rumore equivalente di inngresso al TIA

%Prx_bit1=ptx_bit1 * (A_fspecial)^2;  %potenza ricevuta nel caso di compensazione
%Prx_bit0=ptx_bit0 * (A_fspecial)^2; %potenza ricevuta nel caso di compensazione
Prx_bit1 = ptx_bit1 * (A_f)^2; %potenza ricevuta con att. no f.speciale
Prx_bit0 = ptx_bit0 * (A_f)^2; %potenza ricevuta con att. no f.speciale

PI = (rx_intensity) .* Prx_bit1;  %potena istantanea considerando la potenza ricevuta
I = PI .* responsivity;  % fotocorrente[A]
sigma= NEC * sqrt(Bn);  %deviazione stand.del filtro t.impedenza
varianza= sigma.^2;
%w_noise= wgn(1,(LEN*100)+1,varianza,'linear'); %AWGN con varianza sigma.^2
w_noise = sigma*randn(1,(LEN*100)+1); %wnoise ridefinito con randi

figure(6); %rappresentazione del rumore bianco nel tempo
plot(t_s, w_noise );
title('AWGN','color','blue');



%% filtro gaussiano "TIA" 

sigmagauss = 1.2e10; %dev.standard calcolata in modo da avere f_taglio a 10GHz
A = 1; %altezza nell'origine della gaussiana

gaussian = A*exp(-(f_simmetrico).^2 / sigmagauss^2); %definizione gaussiana nel dominio delle frequenze

gaussian_shifted=fftshift(gaussian); %shift della gaussiana a monobanda

% figure(6);  %filtro gaussiano simmetrico
% plot(f_simmetrico, gaussian);  
% grid on;

% figure(7); %filtro gaussiano shiftato con fftshift
% plot(f_s, gaussian_shifted);
% grid on;

%% Passaggio del segnale per il TIA e relativa aggiunta di rumore
rx_awgn= PI + w_noise ; %somma del mio segnale al rum.eq. TIA in ingresso
RX_awgn = fft(rx_awgn); %passaggio nel dominio delle freq.
RX_tia = RX_awgn .* gaussian_shifted; %segnale post filtro gaussiano TIA
rx_tia= ifft(RX_tia); %segnale riportato nel tempo    

figure(10)  %segnale pre tia con rumore
plot(t_s, rx_awgn)
title('Segnale pre TIA con rumore ','color','blue')

figure(7); %segnale con rumore post filtraggio gaussiano  TIA
plot(t_s, real(rx_tia)); 
title('Segnale post TIA ','color','blue');

%% SNR teorico
prx_mediateenuata = (Prx_bit1 + Prx_bit0)/2;
SNR = prx_mediateenuata/sigma^2;
SNR_dB = 10*log(SNR);
mu1 = ptx_bit1 * (A_f)^2;
mu0 = ptx_bit0 * (A_f)^2;
Q = (mu1 - mu0)/ (2*sigma);
P_be = qfunc(Q)

%% Campionamento 
vettorecampioni = zeros(size(bit));  %vettore di campioni estratti
j=0;
for k= 1 : LEN
    vettorecampioni(k) = real(rx_tia((tsimbolo/dt)/2 + j*(tsimbolo/dt)));
    j=j+1;
end

%% Decisore con soglia posta mean

bit_stima = zeros(size(bit));
soglia2 = real(mean(rx_tia));
for scorri= 1: LEN
    if( vettorecampioni(scorri) > soglia2 )
        bit_stima(scorri) = 1 ;
    end
end

%% BER e EyeDiagram 

numerobiterrati = 0;
for scorri= 1: LEN
    if( bit(scorri) ~=  bit_stima(scorri) )
        numerobiterrati= numerobiterrati+1 ;
    end
end

BER = numerobiterrati / LEN

 
vettoretsimbolo= 0 : dt : tsimbolo;  %Con ciclo For, lento con + di 128bit

% for i= 0 : LEN-1
%      for k= 1:101
%      j=k+(100*i);
%      vettoresupporto(k) = rx_awgn(j);
%      end
%      figure(11);
%      title('EYE-diagram ','color','blue');
%      scatter(vettoretsimbolo,vettoresupporto);
%      hold on;
% end
% clear all;

%% funzione rettangolo

function y=rect(x, intervallo, A)

LEN = 16;
bit_intervalli = 100;
y= zeros(size(x));
    for i = 1:LEN
         if i >= intervallo && i < (intervallo +A)
            y((((i-1)*bit_intervalli)+1) : (i*bit_intervalli +1)) = 1;
         end
    end    


% A ampiezza del rettangolo
% centro

end
