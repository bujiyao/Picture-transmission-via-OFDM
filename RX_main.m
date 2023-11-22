clear all,clc,close all;
warning('off','plutoradio:sysobj:FirmwareIncompatible');
folder = fileparts(which(mfilename)); 
addpath(genpath(folder))
count = 1;
*
H = 120; W = 160;       % 圖片大小         
CH = 1;                 % 黑白=1
                        % 彩色=3

Total_bits = H*W*8*CH;  % 設定全部的資料長度

HS_total = 10;           % 設定資料要拆成幾組做傳送

Part_bits = Total_bits/HS_total;    % 設定每一組的資料長度
RS_intervel = 2;         % 設定RS間隔

N_cpc = 4;               % 設定調變

fft_size = 256;          % 設定子載波大小
CRC_Type = 241;          % CRC的種類
CRC_SIZE = 24;           % CRC的大小
code_rate = 1/3;

Rx_gain = 50;            % 設定rx的功率大小
Tx_gain = -20;           % 設定feedback_tx的功率大小

fc = 3500e6;               % 設定data的中心頻率
fc_HS = 2000e6;            % 設定feedback的中心頻率

fs = 60e6;                 % 設定採樣頻率

test_data = randi([0 1],1,Total_bits);
% fid = load("3_120160.txt","r");  % 讀檔
% test_data = fid;
%%%%%%%%%%%%%%%%%%%%%%%%   Index   %%%%%%%%%%%%%%%%%%%%%%%
Inx_carrier = [1:64,fft_size-64+1:fft_size];    % 設定data佔用子載波的位置
Inx_carrier_sync = [1:63,fft_size-63:fft_size];  % 設定sycn bits佔用子載波的位置

test_data = reshape(test_data,Part_bits,[]);   % 將資料分成n組
test_data = test_data.';
test_data_CRC_num = Part_bits + CRC_SIZE;
test_data_encoding_num = test_data_CRC_num/code_rate;

RS_num = fix(test_data_encoding_num/N_cpc/length(Inx_carrier)/(RS_intervel-1)) + 1 ;  % RS的數量 
Num_Data_Bit = RS_num*N_cpc*length(Inx_carrier)*(RS_intervel-1);
while(mod(Num_Data_Bit,1/code_rate) ~= 0 || mod(Num_Data_Bit-test_data_encoding_num,1/code_rate) ~= 0)
    RS_num = RS_num + 1;
    Num_Data_Bit = RS_num*N_cpc*length(Inx_carrier)*(RS_intervel-1);
end

Data_bit = [test_data, test_data(:,end+1-(Num_Data_Bit-test_data_encoding_num)*code_rate:end)];
RS_end = RS_num*RS_intervel + 2;  % 最後一個symbol的位置

Inx_RS = 2:RS_intervel:RS_end;    % RS從第二個symbol
Inx_data = 3 :(Inx_RS(1,end)-1);  % data從第三個symbol
for i =1:length(Inx_RS)-1
    Inx_temp = find(Inx_data == Inx_RS(i));
    Inx_data(Inx_temp) = [];
end

Inx_RS = 2:RS_intervel:RS_end;    % RS從第二個symbol
Inx_data = 3 :(Inx_RS(1,end)-1);  % data從第三個symbol
for i =1:length(Inx_RS)-1
    Inx_temp = find(Inx_data == Inx_RS(i));
    Inx_data(Inx_temp) = [];
end

Total_len = length(Inx_data) + length(Inx_RS) + 1;
CP_length = 32;        % CP佔用的子載波數量

tx_data_freq = zeros(fft_size,Total_len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%   PSS   %%%%%%%%%%%%%%%%%%%%%%%%%
N1 = 0;
N2 = 0;
[CellID,PSS,SSS] = SS(N1,N2);
tx_data_freq(Inx_carrier_sync,1) = PSS; % 將同步位元放在第一個symbol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   CSIRS   %%%%%%%%%%%%%%%%%%%%%%%%
for k=1:numel(Inx_RS)
    c= Inx_RS(1,k);
    n=mod(floor(c/14)+1,10);  % is the slot number within a radio frame
    l=mod(c,14);              % is the OFDM symbol within a slot
    RS_Seq = CSIRS(n,l,length(Inx_carrier));
    tx_data_freq(Inx_carrier, Inx_RS(1, k)) = RS_Seq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   data   %%%%%%%%%%%%%%%%%%%%%%%%%
final_data = zeros(size(test_data));

filename_total = "result"; %總資料夾
if exist(filename_total,"dir")
    rmdir(filename_total, 's')
end
mkdir(filename_total)
while true
    % 建立常用資料夾
    filename_count = string(count);
    filename_count_fr = filename_total + '\' + filename_count + '\fr';
    mkdir(filename_count_fr)
    filename_count_dp = filename_total + '\' + filename_count + '\dp';
    mkdir(filename_count_dp)
    filename_count_constellation = filename_total + '\' + filename_count + '\constellation';
    mkdir(filename_count_constellation)
    filename_count_camera = filename_total + '\' + filename_count + '\camera';
    mkdir(filename_count_camera)

    for HS_index = 1:HS_total
        close all;
        filename = string(HS_index) + '.png';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%      CRC       %%%%%%%%%%%%%%%%%%%
        Data_Bit = double(CRC_Gen_NR_V2(Data_bit(HS_index,:), CRC_Type)); % add CRC bits
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%   turbo code   %%%%%%%%%%%%%%%%%%%
        % Variables
        info_len = length(Data_Bit); % Frame length
        max_iteration = 18; % Turbo decoder iteration upper bound
        Ec_No_dB = 100000; % Channel SNR
    
        % Generate RSC encoders
        transitions = polynomial2trellis([[1 1] ; [5 7]]);
    
        % Generate information sequence
        info_seq = Data_Bit; 
        
        % Generate interleaver and de-interleaver
        block = fix(sqrt(info_len)) + 1;
        interleaver = [1:info_len -ones(1,block^2-info_len)];
        interleaver = reshape(rot90(reshape(interleaver,block,block)',2),1,[]);
        temp_index = find(interleaver == -1);
        interleaver(temp_index) = [];
    
        % Generate select matrix
        select_matrix = [[1 1] ; [1 1]];
            
        % Conduct channel coding
        encoded_seq = turbo_encoder(info_seq, transitions, interleaver, select_matrix);
    %     encoded_seq = info_seq;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%   modulation   %%%%%%%%%%%%%%%%%%%
        tx_data_mod = modulation_LTE(encoded_seq,N_cpc);
        tx_data_mod = reshape(tx_data_mod,numel(Inx_carrier),[]);
        tx_data_freq(Inx_carrier,Inx_data) = tx_data_mod;
        tx_data = ifft(tx_data_freq,[],1)*sqrt(fft_size);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%   CP   %%%%%%%%%%%%%%%%%%%%%%%%%%%
        tx_data_CP = zeros(fft_size + CP_length,Total_len);
        tx_data_CP((1 : CP_length), :) = tx_data((end - CP_length + 1 : end), :);
        tx_data_CP((CP_length +1 : end), :) = tx_data(:, :);
        sync_time = tx_data_CP(:,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%   plutoSDR   %%%%%%%%%%%%%%%%%%%%%%%%%%
        Ts = 1/fs;
        tx_data_seq = reshape(tx_data_CP,[],1);
        % [rx_radio,tx_radio]=set_pluto(fc,fs,sn,Rx_Gain,Tx_Gain)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%   CRC 重傳   %%%%%%%%%%%%%%%%%%%%%%%%%%
        CRC_result = 0;
        while (CRC_result ==  0)
            close all;
            rx_radio=set_pluto_RX(fc,fs,tx_data_seq,Rx_gain);
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%   時間同步   %%%%%%%%%%%%%%%%%%%%%%%%%%%
            while true
                [rx_data_seq, valid, overflow] =  rx_radio();
        
                [a, lags] = xcorr(rx_data_seq,sync_time); %交相關
                
    %             figure('Name','時間同步:corr抓時間','NumberTitle','off');
    %             subplot(2,1,1);
    %             plot(lags,abs(a));grid on; %tx rx交相關圖
                
                %利用最高值，找出每個符合的時間
                XX = max(abs(a))*0.4;
                [TF2,P2] = islocalmax(abs(a),'MinProminence',XX);
                
                %用來判斷是否傳了下一組
                if(length(lags(TF2)) == 2) 
                    break;
                end
            end
        
            %釋放tx rx接口
            release(rx_radio);
        
            %找出最小的td，對rx_data做時間同步
            all_td = lags(TF2);
            b = find(all_td > 0);
            td = all_td(min(b));
            rx_data_seq = rx_data_seq(td+1:end - (length(tx_data_seq)-td) );
            
            %確認是否成功同步，其td為0
            [a , lags] = xcorr(rx_data_seq,sync_time);
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
            
                %將另一個RS的變數存起來
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
            
            % 頻域通道
            f = linspace(-fs/2,fs/2,fft_size);
            figure('Name','頻域:channel','NumberTitle','off');
            
            plot(f,abs(fftshift(CE_freq(:,2))));
            xlabel("Sampling Freq.")
            ylabel("Amplitude")
            xticks([-fs/2 -fs/4 0 fs/4 fs/2]);
            xticklabels({'fc-BW/2','fc-BW/4','fc','fc+BW/4','fc+BW/2'});
            ylim([0 0.2]);
            filename_fr = filename_total + '\' + filename_count + '\fr\' + filename;
            saveas(gcf,filename_fr)

            % delay profile
            t = linspace(0,1/fs,fft_size);
            h=ifft(CE_freq);
            figure('Name','時域:channel','NumberTitle','off');
            plot(t,abs(h(:,2)))
            xlabel("Time")
            ylabel("Power")
            ylim([0 0.025]);
            filename_dp = filename_total + '\' + filename_count + '\dp\' + filename;
            saveas(gcf,filename_dp)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%   等化器後   %%%%%%%%%%%%%%%%%%%%%%%%%
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

            filename_constellation = filename_total + '\' + filename_count + '\constellation\' + filename;
            saveas(gcf,filename_constellation)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%  remove turbo %%%%%%%%%%%%%%%%%%%%%%%
            decoded_seq = double(turbo_decoder(infor_freq, transitions, interleaver, ...
                select_matrix, max_iteration, Ec_No_dB));
    %         decoded_seq = infor_freq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%   CRC + BER   %%%%%%%%%%%%%%%%%%%%%%%
            fprintf("\n第 %i 組BER",HS_index);
            CRC_result = CRC_Check_NR_V2(decoded_seq,CRC_Type) % check CRC
            decoded_seq = decoded_seq(1:end-24);              % remove CRC bits
            final_data_1 = decoded_seq(1:length(test_data));   % remove 尾端資料
            BER = biterr(final_data_1,test_data(HS_index,:))/length(test_data(HS_index,:))
            final_data(HS_index,:) = final_data_1;
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%   Feedback   %%%%%%%%%%%%%%%%%%%%%%%%
        if(mod(HS_index,2) == 1)
            feedback_signal = complex(PSS);     %feedback傳的資料
            tx_data_freq(Inx_carrier_sync,1) = SSS;      %下組資料用的同步訊號，用來判斷是否已傳下一組
        else
            feedback_signal = complex(SSS);
            tx_data_freq(Inx_carrier_sync,1) = PSS;
        end
    
        %%%%%%%%%%%%%%%%%%%%   一台兩台差別   %%%%%%%%%%%%%%%%%%%%%
        [tx_radio_HS]=set_pluto_TX_HS(fc_HS,fs,feedback_signal,Tx_gain);
        tx_radio_HS.transmitRepeat(feedback_signal.');
    %     pause;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        release(rx_radio);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%   BER   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BER_total = biterr(final_data,test_data)/numel(test_data)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%   picture   %%%%%%%%%%%%%%%%%%%%%%%%%
    final_data_seq = reshape(final_data.',[],1);
    figure('Name','圖片','NumberTitle','off')
    bit_2_img(final_data_seq.',H,W,CH);
    
    filename_camera = filename_total + '\' + filename_count + '\camera\' + filename;
    saveas(gcf,filename_camera)
    filename_camera = filename_total + '\' + filename_count + '\all_data.mat';
    save(filename_camera);

    count = count + 1;

end