%%-- 2º Trabalho Computacional - Análise de Sinais Biológicos ------------%
%--- Nome: Amanda Rosa F. Jorge - RRA: 225316 ----------------------------%

clc;
clear all;
close all;

%% -------------------------- Exercicio 1 --------------------------------%

%----------------------- Seleção de Neurônios ----------------------------%

matObj = matfile('tc2ex1.mat');
Documento = fopen('Resultados.txt','w');

spikes = matObj.spikes;
time = matObj.time;
force = matObj.force;

%Configurar matriz tempo-força
    
    Matriz_forca_tempo(:,1) = time(:,1);
    Matriz_forca_tempo(:,2) = force(:,1);

vetor_neuron = zeros(1,10);

fprintf(Documento, 'Neuronios Selecionados\n\n');
fprintf(Documento, 'Coluna 1 - Índices Neurônios e Coluna 2 - Quantidade de Spikes\n\n');

for i = 1:10 

    neuron_rand = randi([1 100],1);
    
    index_neuron = min(find(neuron_rand == spikes(:,2)));
    apontador = spikes(index_neuron,2);
    contador = 1;
    
    while(spikes(index_neuron,2) == neuron_rand)
        
        contador = contador + 1;
        index_neuron = index_neuron + 1;
        
    end
    
    if (contador>= 50)
        
        if(isempty(find(neuron_rand == vetor_neuron)) ^ neuron_rand ~= 0 )
        
        vetor_neuron(1,i) = neuron_rand; %coleta índice do neurônio
        vetor_neuron(2,i) = contador-1; %coleta quantidade de disparos do neurônio selecionado
        fprintf(Documento,'%d -> ',vetor_neuron(1,i));
        fprintf(Documento,'%d\n',vetor_neuron(2,i));
        
        end
    end
end

% %------------------- NEURÔNIOS SELECIONADOS ----------------------------%
% 
% %--------------------------- Força média -------------------------------%

%removendo o transitório de forca de 1 a 7000 amostras
Force = force(7000:end);
Media_forca = mean(Force);
STD_forca = std(Force);
CV_forca = STD_forca/Media_forca;

fprintf(Documento, '------------------------------------------------\n\n');
fprintf(Documento, 'Média de Força(N): %f\n', Media_forca);
fprintf(Documento, 'Desvio Padrão de Força: %f\n', STD_forca);
fprintf(Documento, 'Coeficiente de Variação de Força: %f\n', CV_forca);

% Matriz Instantes Spikes
 for n = 1:10   
    neuron = vetor_neuron(1,n);
    apontador = min(find(neuron == spikes(:,2)));
    disparos = vetor_neuron(2,n);
    
    for instante = 1:disparos
        
        Matriz_time(instante,n) = spikes(apontador,1);
        
        if(spikes(apontador,1)>=350) %desconsiderando os primeiros 350ms de transitório de força
            Matriz_time_est(instante,n) = spikes(apontador,1); %matriz tempos linha tempos e coluna neuronios
            apontador = apontador + 1;
        else
            Matriz_time_est(instante,n) = 0;
            apontador = apontador + 1;
        end
    end
 end
 
 %construindo celulas com os instantes spikes selecionados sem o transitório
 
 celula = cell(10,1);

for n = 1:10
    
    auxiliar = nonzeros(Matriz_time_est(:,n));
    auxiliar = auxiliar';
    celula_est(n,1) = {auxiliar}; 
    
    auxiliar2 = nonzeros(Matriz_time(:,n));
    auxiliar2 = auxiliar2';
    celula(n,1) = {auxiliar2};
end

% %---------------------Intervalos Inter-Spikes---------------------------%
fprintf(Documento, '------------------------------------------------\n\n');
fprintf(Documento, '---------Intervalos Inter-spikes----------------\n\n');

 for n = 1:10
     
     fprintf(Documento, '------------\n');
     fprintf(Documento, 'Neurônio: %i\n\n',n);

     Vetor_tempo = cell2mat(celula(n,1));
     disparos = length(Vetor_tempo);
     
     Vetor_tempo_est = cell2mat(celula_est(n,1));
     disparos_est = length(Vetor_tempo_est);

     ISI = diff(Vetor_tempo);
     
     ISI_medio(1,n) = mean(ISI);
     fprintf(Documento, 'ISI_medio: %d\n\n',ISI_medio(1,n));
     
     %%% os desvios-padrão e os coeficientes de variação dos ISIs para os 10 neurônios selecionados
     DP_ISI(1,n) = std(ISI);
     fprintf(Documento, 'Desvio Padrao ISI: %d\n\n',DP_ISI(1,n));
     
     CV_ISI(1,n) = DP_ISI(1,n)/ISI_medio(1,n); 
     fprintf(Documento, 'Coeficiente de Variacao ISI: %d\n\n',CV_ISI(1,n));
     %%% os coeficientes de assimetria e curtose dos ISIs para os 10 neurônios selecionados
     Assimetria(1,n) = skewness(ISI);
     fprintf(Documento, 'Assimetria ISI: %d\n\n',Assimetria(1,n));
     
     Curtose(1,n) = kurtosis(ISI); 
     fprintf(Documento, 'Curtose ISI: %d\n\n',Curtose(1,n));
      
     %%% os histogramas dos ISIs para os 10 neurônios selecionados
     norm = 100;
     subplot(2,5,n);
     Histograma(1,n) = histogram(ISI); 
     xlim([0 300]);
     title(['Neurônio',num2str(n)]);
     
 end

% % --------------------------Vetor Spike Train---------------------------%
Fs = (length(time)/(time(length(time),1)/1000));
Ts = 1/Fs;
Passo_int = 1/Fs;

% frequencia instantena de spikes considerando o transitorio
for n = 1:10
    
    disparos = vetor_neuron(2,n);
    Vetor_MUAPT(n,1:200000) = zeros(1,length(Matriz_forca_tempo(:,1)));
    index_spike = 0;
    
    for j = 1:disparos
    
    apontador = Matriz_time(j,n);
    index_spike = find(apontador == Matriz_forca_tempo);
%     Vetor_MUAPT(n,index_spike) = 1/Passo_int;
    Vetor_MUAPT(n,index_spike) = 1;

    apontador_est = Matriz_time_est(j,n);
    index_spike_est = find(apontador_est == Matriz_forca_tempo);
    Vetor_MUAPT_est(n,index_spike_est) = 1/Passo_int;
    
    end
    
    Spike_train(n,:) =  Vetor_MUAPT(n,1:200000);
    Spike_train_est(n,:) = Vetor_MUAPT_est(n,1:200000);
    
end

figure();
subplot(10,1,1);
plot(time,Spike_train(1,:));
hold on
title('Vetor Spike Train');
xlim([0 10000]);

for n = 2:10
    
    subplot(10,1,n);
    plot(time,Spike_train(n,:));
    hold on
end

xlabel('Tempo(ms)');
ylabel('Disparo');

% % ----------------------- Estimação Frequência Instantânea --------------%
fprintf(Documento, '------------------------------------------------\n\n');
fprintf(Documento, '---------Estimativa Frequência Instantânea------\n\n');

L = length(time)*(0.05); % janela de 500 ms
w = hann(L);
area = cumtrapz(w);
janela_unitaria = w./area(end,1);
    
%considerando o transitorio de força
for p = 1:10
    
    fprintf(Documento, '------------\n');
    fprintf(Documento, 'Neurônio: %i\n\n',p); 
    Freq_inst(:,p) = conv(Spike_train(p,:),janela_unitaria);
    Freq_inst_media(p,1) = mean(Freq_inst(10000:190000,p)); 
    fprintf(Documento, 'Frequencia Instantânea Média (Hz): %f\n\n',Freq_inst_media(p,1));
    
end  

% % ---------------- Analise Estatística ---------------------------------%

%utilizando trechos estacionários, desconsiderando o transitório de força
for n = 1:10
    
    Freq_inst_est(:,n) = conv(Spike_train_est(n,:),janela_unitaria, 'same');
    Freq_inst_media_est(n,1) = mean(Freq_inst_est(10000:180000,n));
    
end   

%removendo a tendência linear
Forca_novo = detrend(Force); %remoção de tendências lineares do sinal força

% Funções correlação cruzada entre frequencias estimadas e sinal força
i=1;
for n = 1:10
    figure();
    aux1_a = Freq_inst_est(10000:180000,n);
    aux2_a = detrend(aux1_a);
    [c_a(i,:),lags_a(i,:)] = xcorr(aux2_a,Forca_novo(10000:end-13001,:),[], 'coeff');
    lags_a(i,:) = lags_a(i,:)./Fs;
    plot(lags_a(i,:),c_a(i,:));
    title([' Funcao Correlacao Cruzada entre Neuronio ', num2str(n),'  e Força']);
    ylabel('Função de Correlacao Cruzada');
    xlabel('Time Lags (s)');
    i=i+1;
    
end

CC_a_mean = mean(c_a);
figure();
plot(lags_a,CC_a_mean); 
title([' Funcao Correlacao Cruzada entre Frequencia Instantânea dos Neuronios e Força']);
ylabel('Função de Correlacao Cruzada');
xlabel('Time Lags (s)');

% Funções correlação cruzada entre pares de frequências estimadas
i=1;
for n = 1:9
    
    for j = n+1:10
        
        aux1_1 = Freq_inst_est(10000:180000,n) - Freq_inst_media_est(n,1);
        aux1_2 = detrend(aux1_1);

        aux1_2 = Freq_inst_est(10000:180000,j) - Freq_inst_media_est(j,1);
        aux2_2 = detrend(aux1_2);

        [c_b(i,:),lags_b(i,:)] = xcorr(aux1_2,aux2_2,[], 'coeff');
        lags_b(i,:) = lags_b(i,:)./Fs;
        i=i+1;
        
    end
    
end

i=1;
for n = 1:9
    figure();
    for j = n+1:10
       if(n==1) %10 graficos
           place = i;
            subplot(5,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==2) %9 graficos
            place = i-9;
            subplot(5,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==3) %8 gráficos
             place = i-17;
            subplot(4,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==4) %7 gráficos
            place = i-24;
            subplot(4,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==5)%6 gráficos
             place = i-30;
            subplot(3,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==6)%5 gráficos
             place = i-35;
            subplot(3,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==7)%4 gráficos
            place = i-39;
            subplot(2,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==8)%3 gráficos
            place = i-42;
            subplot(2,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        elseif(n==8)%2 gráficos
            place = i - 44;
            subplot(1,2,place);
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
        else %1 gráfico
            plot(lags_b(i,:),c_b(i,:));
            title([' Funcao Correlacao Cruzada entre Pares de Neuronios ', num2str(n), ' e ', num2str(j)] );
            ylabel('Função Corr. Cruzada');
            xlabel('Time Lags (s)');
            
       end
        i=i+1;
    end
end

fclose(Documento);
