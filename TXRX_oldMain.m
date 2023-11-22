clear all,clc,close all
warning('off','plutoradio:sysobj:FirmwareIncompatible');
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

Total_bits = 120*160*8*3;
RS_intervel = 3;

N_cpc = 2;
fft_size = 256;

fc=1200e6;
fs=60e6;
%%%%%%%%%%%%%%%%%%%%%%%%   Index   %%%%%%%%%%%%%%%%%%%%%%%
Inx_carrier = [1:80,fft_size-80+1:fft_size];
Inx_carrier_PSS = [1:63,fft_size-63:fft_size];

RS_end = Total_bits/N_cpc/length(Inx_carrier)*RS_intervel/(RS_intervel-1)+2;

Inx_RS = 2:RS_intervel:RS_end; 
Inx_data = 3 :(Inx_RS(1,end)-1);
for i =1:length(Inx_RS)-1
    Inx_temp = find(Inx_data == Inx_RS(i));
    Inx_data(Inx_temp) = [];
end

Total_len = length(Inx_data) + length(Inx_RS) + 1;
CP_length =32;

tx_data_freq = zeros(fft_size,Total_len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   PSS   %%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 0;
N2 = 0;
[CellID,PSS,SSS] = SS(N1,N2);
tx_data_freq(Inx_carrier_PSS,1) = PSS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   CSIRS   %%%%%%%%%%%%%%%%%%%%%%%%
for k=1:numel(Inx_RS)
    c= Inx_RS(1,k);
    n=mod(floor(c/14)+1,10);  %is the slot number within a radio frame
    l=mod(c,14);              %is the OFDM symbol within a slot
    RS_Seq = CSIRS(n,l,length(Inx_carrier));
    %RS_Seq_mod = modulation_LTE(RS_Seq,N_cpc);
    tx_data_freq(Inx_carrier, Inx_RS(1, k)) = RS_Seq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   data   %%%%%%%%%%%%%%%%%%%%%%%%%
% test_data = randi([0 1],1,numel(Inx_carrier)*length(Inx_data)*N_cpc);
% fid = load("picture\0606_BIN_DATA_GRAY.txt","r");
fid = load("3_120160.txt","r");  % 讀檔
test_data = fid;

tx_data_mod = modulation_LTE(test_data,N_cpc);
tx_data_mod = reshape(tx_data_mod,numel(Inx_carrier),[]);
tx_data_freq(Inx_carrier,Inx_data) = tx_data_mod;
tx_data = ifft(tx_data_freq,[],1)*sqrt(fft_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   CP   %%%%%%%%%%%%%%%%%%%%%%%%%%%
tx_data_CP = zeros(fft_size + CP_length,Total_len);
tx_data_CP((1 : CP_length), :) = tx_data((end - CP_length + 1 : end), :);
tx_data_CP((CP_length +1 : end), :) = tx_data(:, :);
PSS_time = tx_data_CP(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   plutoSDR   %%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 1/fs;
tx_data_seq = reshape(tx_data_CP,[],1);
% tx_data_seq = reshape(tx_data,[],1);
% [rx_radio,tx_radio]=set_pluto(fc,fs,sn,Rx_Gain,Tx_Gain)
[rx_radio,tx_radio]=set_pluto(fc,fs,tx_data_seq,30,-50);
tx_radio.transmitRepeat(tx_data_seq);
[rx_data_seq, valid, overflow] =  rx_radio();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   時間同步   %%%%%%%%%%%%%%%%%%%%%%%%%%%
[a, lags] = xcorr(rx_data_seq,PSS_time); %交相關

figure('Name','時間同步:corr抓時間','NumberTitle','off');
subplot(2,1,1);
plot(lags,abs(a));grid on; %tx rx交相關圖

%利用最高值，找出每個符合的時間
XX = max(abs(a))*0.5;
[TF2,P2] = islocalmax(abs(a),'MinProminence',XX);

subplot(2,1,2);
stem(lags(TF2),abs(a(TF2))); %符合時間點圖

%找出最小的td，對rx_data做時間同步
all_td = lags(TF2);
b = find(all_td > 0);
td = all_td(min(b));
rx_data_seq = rx_data_seq(td+1:end - (length(tx_data_seq)-td) );


%確認是否成功同步，其td為0
[a , lags] = xcorr(rx_data_seq,PSS_time);
figure('Name','時間同步:確認是否正確同步','NumberTitle','off');
plot(lags,abs(a));grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   移除CP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx_data_cp = reshape(rx_data_seq,fft_size + CP_length,[]);
rx_data = rx_data_cp((CP_length + 1 :end),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   頻率偏移補償   %%%%%%%%%%%%%%%%%%%%%%%
rx_data_freq = fft(rx_data,[],1)/sqrt(fft_size);
rx_data_aft_seq = rx_data_seq;
rx_data_aft_seq_RS = zeros(fft_size+CP_length,size(Inx_RS,2));

for i = 1:(size(Inx_RS,2)-1)
        
    rn_PramblePSS_1 = rx_data(:,Inx_RS(1,i));
    rn_PramblePSS_1 = rn_PramblePSS_1./tx_data(:,Inx_RS(1,i));
        
    rn_PramblePSS_2 = rx_data(:,Inx_RS(1,i+1));
    rn_PramblePSS_2 = rn_PramblePSS_2./tx_data(:,Inx_RS(1,i+1));

    rn_freq_pramble_phasePSS = rn_PramblePSS_1.*conj( rn_PramblePSS_2 ); %找出兩個PSS之相位差
    Freq_OffetPSS=atan2(imag(rn_freq_pramble_phasePSS),real(rn_freq_pramble_phasePSS));
    fe_est_PSS =mean(Freq_OffetPSS)/(2*pi*Ts*(fft_size+CP_length)*(Inx_RS(1,i+1) - Inx_RS(1,i)));

    %對rx_data做頻率補償
    temp_inx = ((fft_size+CP_length)*(Inx_RS(1,i)-1)+1):(fft_size+CP_length)*(Inx_RS(1,i+1));
    rx_data_aft_seq(temp_inx) = rx_data_seq(temp_inx) ...
                       .* exp(1i*2*pi*fe_est_PSS*Ts*(0:length(rx_data_seq(temp_inx)) - 1)).';

%     temp_inx2 = ((fft_size+CP_length)*(Inx_RS(1,i)-1)+1):(fft_size+CP_length)*(Inx_RS(1,i));
%     temp_inx3 = ((fft_size+CP_length)*(Inx_RS(1,i+1)-1)+1):(fft_size+CP_length)*(Inx_RS(1,i+1));
%     rx_data_aft_seq(temp_inx2) = rx_data_aft_seq(temp_inx2)/2 + rx_data_aft_seq(temp_inx3)/2;
    %對RS做平均
    rx_RS_temp_inx2 = ((fft_size+CP_length)*(Inx_RS(1,i+1)-1)+1):(fft_size+CP_length)*(Inx_RS(1,i+1));
    rx_data_aft_seq_RS(:,i) = rx_data_aft_seq(rx_RS_temp_inx2);
end
%移除cp
rx_data_aft_cp = reshape(rx_data_aft_seq,fft_size + CP_length,[]);
rx_data_aft = rx_data_aft_cp((CP_length + 1 :end),:);
rx_data_aft_seq_RS = rx_data_aft_seq_RS((CP_length + 1 :end),:);

%轉到頻域
rx_data_aft = rx_data_aft(:,(1:Total_len));
rx_data_aft_freq = fft(rx_data_aft,[],1)/sqrt(fft_size);
rx_data_aft_seq_RS_freq = fft(rx_data_aft_seq_RS,[],1)/sqrt(fft_size);

%查看星座圖的分布
figure('Name','頻偏:data','NumberTitle','off');
plot( real( rx_data_aft_freq(Inx_carrier,Inx_data(1,1)) ),imag( rx_data_aft_freq(Inx_carrier,Inx_data(1,1)) ),"b*",'MarkerSize',10);
hold on;
plot( real( rx_data_aft_freq(Inx_carrier,Inx_data(1,2)) ),imag( rx_data_aft_freq(Inx_carrier,Inx_data(1,2)) ),"r*",'MarkerSize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   頻譜分析   %%%%%%%%%%%%%%%%%%%%%%%%
%兩個取平均
CE_freq = zeros(fft_size, Total_len);
CE_freq_RS = zeros(fft_size,size(Inx_RS,2));
for i = 1:size(Inx_RS,2)-1
    CE_freq(Inx_carrier,Inx_RS(1,i)) = rx_data_aft_freq(Inx_carrier,Inx_RS(1,i))./tx_data_freq(Inx_carrier,Inx_RS(1,i));
    CE_freq_RS(Inx_carrier,i) = rx_data_aft_seq_RS_freq(Inx_carrier,i)./tx_data_freq(Inx_carrier,Inx_RS(1,i+1));

    CE_avg = CE_freq(Inx_carrier,Inx_RS(1,i))/2 + CE_freq_RS(Inx_carrier,i)/2;

    for n = (Inx_RS(1,i)+1 : Inx_RS(1,i+1)-1)
        CE_freq(Inx_carrier,n )  = CE_avg;
    end
end


% 通道在頻率(Inx_RS = [2])
f = linspace(0,fs,fft_size);
figure('Name','頻域:channel','NumberTitle','off');
plot(f,abs(CE_freq(:,Inx_RS)));

% 通道在時域(Inx_RS = [2])
t = linspace(0,1/fs,fft_size);
h=ifft(CE_freq);
figure('Name','時域:channel','NumberTitle','off');
plot(t,abs(h(:,Inx_RS)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   等化器前   %%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','(前)頻域: information','NumberTitle','off');
rx_data_aft_data = rx_data_aft_freq(Inx_carrier,Inx_data);
rx_data_aft_data_seq = reshape(rx_data_aft_data,[],1);
plot( real( rx_data_aft_data_seq ),imag( rx_data_aft_data_seq ),"b*",'MarkerSize',10);

% rx_data_demod = demodulation_LTE(rx_data_aft_data_seq,N_cpc);
% [index, mapper_binary]  = demap_test(rx_data_demod.',N_cpc);
% infor_freq_before = zeros(length(index),N_cpc);
% for i = 1:numel(index)
%     infor_freq_before(i,:) =mapper_binary(index(i,:),:) ;
%     infor_freq_before = [infor_freq_before;];
% end
% infor_freq_before = reshape(infor_freq_before', 1, numel(infor_freq_before));
% 
% BER = biterr(infor_freq_before,test_data)/length(test_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   等化器後   %%%%%%%%%%%%%%%%%%%%%%%%
rx_equal = rx_data_aft_freq(Inx_carrier,Inx_data)./CE_freq(Inx_carrier,Inx_data);
rx_equal_seq = reshape(rx_equal,[],1);
figure('Name','頻域: information','NumberTitle','off');
plot( real( rx_equal_seq ),imag( rx_equal_seq ),"b*",'MarkerSize',10);
axis([-2 2 -2 2]);

rx_data_demod = demodulation_LTE(rx_equal_seq,N_cpc);

[index, mapper_binary]  = demap_test(rx_data_demod.',N_cpc);
infor_freq = zeros(length(index),N_cpc);
for i = 1:numel(index)
    infor_freq(i,:) =mapper_binary(index(i,:),:) ;
    infor_freq = [infor_freq;];
end
infor_freq = reshape(infor_freq', 1, numel(infor_freq));

BER = biterr(infor_freq,test_data)/length(test_data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','圖片','NumberTitle','off')
bit_2_img(infor_freq,160,120,3);

