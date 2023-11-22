function [index,mapper_binary ]= demap_test(rn_demod,N_cpc_Data)   % insert the demodulated I/Q data
% load('mapper_IQ_256QAM');                   
% load('mapper_binary_256QAM');       

k_BPSK=1;
k_QPSK=1;
k_16QAM=1;
k_64QAM=1;
k_256QAM=1;
k_1024QAM=1;
switch(N_cpc_Data)

    case 1 % BPSK
    mapper_IQ =cell2mat(struct2cell(load('mapper_IQ_BPSK.mat')));                 
    mapper_binary = cell2mat(struct2cell(load('mapper_binary_BPSK.mat')));  

        
    case 2 % QPSK
    mapper_IQ   = cell2mat(struct2cell(load('mapper_IQ_QPSK.mat')));                 
    mapper_binary = cell2mat(struct2cell(load('mapper_binary_QPSK.mat')));  

        
    case 4 % 16QAM
    mapper_IQ  =   cell2mat(struct2cell(load('mapper_IQ_16QAM.mat')));               
    mapper_binary = cell2mat(struct2cell(load('mapper_binary_16QAM.mat')));


        
    case 6 % 64QAM
    mapper_IQ =cell2mat(struct2cell(load('mapper_IQ_64QAM.mat')));                  
    mapper_binary=cell2mat(struct2cell(load('mapper_binary_64QAM.mat')));

        
    case 8 % 256QAM
    mapper_IQ =  cell2mat(struct2cell(load('mapper_IQ_256QAM.mat'))).';                  
     mapper_binary =    cell2mat(struct2cell(load('mapper_binary_256QAM.mat')));

    case 10 % 1024QAM
    mapper_IQ   = cell2mat(struct2cell(load('mapper_IQ_1024QAM.mat'))).';             
    mapper_binary = cell2mat(struct2cell(load('mapper_binary_1024QAM.mat'))); 
           
end

index= rn_demod == mapper_IQ;
[index,~]=find(index);
end