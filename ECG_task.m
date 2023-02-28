%% Trabalho 1 - Disciplina Análise de Sinais Biológicos
% Prof: Dr. Leonardo Abdala
% Nome: Amanda Rosa Ferreira Jorge


%% Sinal Original
file = 'C:\Users\amand\OneDrive\Documentos\MATLAB\Analise de sinais\1 Trabalho\signal.txt';

signal = fopen(file,'r');

A = fscanf(signal,'%f');

%Definições do sinal
Fs = 500; %frequencia de amostragem
T = 1/Fs; %periodo
L = 10000; %tamanho do sinal em amostras
t = (0:L-1)*T;  % Time vector

%Gráfico do sinal ECG original
figure();
plot(t,A,'-');
legend('ECG signal Original');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');
title(' sinal ECG Original');

%Fast Fourier Transform do sinal Original
FFT = fft(A);

%Espectro calculado em 2 lados (P2) e avaliados de acordo com o tamanho L do
%sinal
P2 = abs(FFT/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%Gráfico da Densidade Espectral do Sinal
%Definindo a frequência e plotando a amplitude do espectro em um lado P1
f = Fs*(0:(L/2))/L;
figure(2);
plot(f,P1);
xlabel ('Frequência(Hz)');
ylabel ('Densidade Espectral de Potência (PSD)');
title(' Análise Espectral do Sinal Original');


%% Filtro Passa-Baixas

wp = 45/250;
ws = 140/250; 
rp = 3;
rs = 40;
[n,Wn] = buttord(wp,ws,rp,rs);
[b,a] = butter(n,Wn, 'low');
dataOut = filter(b,a,A);

%Gráfico do sinal ECG Filtrado Passa-Baixas
figure(3);
plot(t,dataOut,'-');
legend('ECG Filtrado Passa-Baixas');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');
title('Sinal ECG Filtrado Passa-Baixas');

%Fast Fourier Transform do sinal após PB
FFT_1 = fft(dataOut);

%Espectro calculado em 2 lados (P2_1) e avaliados de acordo com o tamanho L do
%sinal
P2_1 = abs(FFT_1/L);
P1_1 = P2_1(1:L/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);

%Gráfico da Densidade Espectral do Sinal
%Definindo a frequência e plotando a amplitude do espectro em um lado P1_1
f = Fs*(0:(L/2))/L;
figure(4);
plot(f,P1_1);
xlabel ('Frequência(Hz)');
ylabel ('Densidade Espectral de Potência (PSD)');

% Função de Transferência Filtro Passa-Baixas
[z,p,k] = butter(n,Wn);
[b,a] = butter(n,Wn);
Hc_s_lowpass = zpk(z,p,k);
[bd,ad] = bilinear(b,a,500);
Hz_lowpass = filt(bd,ad,1/500);

% Resposta em Frequencia Amplitude e Fase do Filtro Passa- Baixas
freqz(b,a,5000,500);
title('Resposta em Frequência do Filtro Passa-Baixas');


%% Filtro Passa-Altas

wp_h = 0.55/250;
ws_h = 0.1/250; 
rp_h = 3;
rs_h = 40;
[n_high,Wn_high] = buttord(wp_h,ws_h,rp_h,rs_h);
[b_high,a_high] = butter(n_high,Wn_high, 'high'); 
dataOut2 = filter(b_high,a_high,dataOut);

%Gráfico do sinal ECG Filtrado Passa-Baixas e Passa-Altas
figure(5);
plot(t,dataOut2,'-');
legend('ECG Filtrado Passa-Baixas e Passa-Altas');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');
title('Sinal ECG Filtrado Passa-Baixas e Passa-Altas');
ylim([-2 2]);
xlim([0 10]);

%Fast Fourier Transform do sinal após PB e PA
FFT_2 = fft(dataOut2);

%Espectro calculado em 2 lados (P2_2) e avaliados de acordo com o tamanho L do
%sinal
P2_2 = abs(FFT_2/L);
P1_2 = P2_2(1:L/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

% Gráfico da Densidade Espectral do Sinal
%Definindo a frequência e plotando a amplitude do espectro em um lado P1_2
f = Fs*(0:(L/2))/L;
figure(6);
plot(f,P1_2);
xlabel ('Frequência(Hz)');
ylabel ('Densidade Espectral de Potência (PSD)');

%Função de Transferência Filtro Passa-Altas
[zh,ph,kh] = butter(n_high,Wn_high,'high');
[bh,ah] = butter(n,Wn);
Hc_s_high = zpk(zh,ph,kh);
[bhd,ahd] = bilinear(bh,ah,500);
Hz_high = filt(bhd,ahd,1/500);

% Resposta em Frequência e Fase do Filtro Passa-Altas
freqz(b_high, a_high,5000,500);
title('Resposta em Frequência do Filtro Passa-Altas');

%% Filtro Notch

% A partir das contas realizadas à mão, encontrou-se a Função de
% Transferência H(z) do Filtro Notch 
a_notch = [0.382];
b_notch = [1 -1.618 1];
dataOut3 = filter(b_notch,a_notch,dataOut2);

%Gráfico do sinal ECG Filtrado Passa-Baixas e Passa-Altas
figure(7);
plot(t,dataOut3,'-');
legend('ECG Filtrado Final');
xlabel('Tempo (s)');
ylabel('Amplitude (mV)');
title('Sinal ECG Filtrado Final');
ylim([-2 2]);
xlim([0 10]);

%Fast Fourier Transform do sinal final
FFT_3 = fft(dataOut3);

%Espectro calculado em 2 lados (P2_3) e avaliados de acordo com o tamanho L do
%sinal
P2_3 = abs(FFT_3/L);
P1_3 = P2_3(1:L/2+1);
P1_3(2:end-1) = 2*P1_3(2:end-1);

% Gráfico da Densidade Espectral do Sinal
%Definindo a frequência e plotando a amplitude do espectro em um lado P1_3
f = Fs*(0:(L/2))/L;
figure(8);
plot(f,P1_3);
xlabel ('Frequncia(Hz)');
ylabel ('Densidade Espectral de Potencia (PSD)');

% Função de Transferência Filtro Notch
Hz_notch = filt(b_notch, a_notch,1/500);

% Resposta em Frequencia Amplitude e Fase do Filtro Notch
freqz(b_notch,a_notch,5000,500);
title('Resposta em Frequência do Filtro Notch');

% %% Gráficos Comparação dos Sinais
% plot(t,A,'-');
% hold on
% plot(t,dataOut,'-');
% hold on
% plot(t,dataOut2,'-');
% legend('ECG Original');
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% ylim([-1.5 1.5]);
% title('ECG em três estados - Original, após PB e após PB e PA');
% xlim([0 10]);
% % %------------------------------%
% plot(f,P1);
% hold on
% % plot(f,P1_1);
% % hold on
% plot(f,P1_3);
% legend('PSD Original');
% xlabel ('Frequência(Hz)');
% ylabel ('Densidade Espectral de Potência (PSD)');
% title('Comparação PSD Original e PSD Final');
% % % % ----------------------------%
% plot(t,A,'-');
% hold on
% plot(t,dataOut,'-');
% hold on
% plot(t,dataOut3,'-');
% legend('ECG Filtrado Final');
% xlabel('Tempo (s)');
% ylabel('Amplitude (mV)');
% ylim([-1.5, 1.5]);
% xlim([0 5]);
% title('Comparação sinal Original e sinal Final');
% 
% % % % % %-------------------------------------------------%
% plot(f,P1);
% hold on
% plot(f,P1_a);
% hold on
% plot(f,P1_3);
% legend('PSD Original');
% xlabel ('Frequência(Hz)');
% ylabel ('Densidade Espectral de Potência (PSD)');
% % % %------------------------------------------------------%
% plot(t,A,'-');
% hold on
% plot(t,dataOut2,'-');
% hold on
 
% 
