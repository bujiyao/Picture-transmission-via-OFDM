%=========================================================================
% 5G NR CRC Generator
% Input: 
%          data: input binary data (row vector)
%          CRC_case: 11 : NR CRC11
%                    16 : NR CRC11
%                    241: NR CRC24A
%                    242: NR CRC24B
%                    243: NR CRC24C
%                    249: LTE CRC
%                    99 : Ethernet
% Output:
%          output: data bit appended with CRC bits
%                              
% Note:  
%       1. The function is implemented by uint32 
%
%   Reference: 
%   [1] 3GPP TS 38.212 V1.1.0 (2017-11), section 5.1
%   Data: 2017.12.25 
%   Author: Zanyu, Chen   
%==========================================================================


function output = CRC_Gen_NR_V2(data, CRC_case)

    switch CRC_case
        case 11  %CRC11, 0xC4200000
            %poly_bin = [1 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            CRC_size = 11;
        case 16  %CRC16, 0x10210000
            %poly_bin = [0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            CRC_size = 16;
        case 241 %CRC24A, 0x864CFB00
            %poly_bin = [1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0];
            CRC_size = 24;
        case 242 %CRC24B, 0x80006300
            %poly_bin = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0 0 0 0 0 0 0];
            CRC_size = 24;
        case 243 %CRC24C, 0xB2B11700
            %poly_bin = [1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1 0 0 0 0 0 0 0 0];
            CRC_size = 24;
        case 249 % LTE_CRC, 0x864CFB00
            %poly_bin = [1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0];
            CRC_size = 24;            
        case 99 %Ethernet: 0x04C11DB7
            %poly_bin = [0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];
            CRC_size = 32;
    end

    tmp = mod(length(data),8);
    CRC_tmp= Compute_CRC32_C_v2(data, CRC_case);
    CRC_tmp = typecast(CRC_tmp, 'uint32');
    CRC = de2bi(CRC_tmp(1), 32, 'left-msb');
    output = [data CRC(1:CRC_size)];        
end

