%=========================================================================
% 5G NR CRC Check
% Input: 
%          data_CRC: input data with CRC (row vector)
%          CRC_case: 11 : NR CRC11
%                    16 : NR CRC11
%                    241: NR CRC24A
%                    242: NR CRC24B
%                    243: NR CRC24C
%                    249: LTE CRC
%                    99 : Ethernet
% Output:
%          CRC_result: 
%              1: CRC check successful
%              0: CRC check failed 
%                              
% Note:  
%       1. The function is implemented by uint32 
%
%   Reference: 
%   [1] 3GPP TS 38.212 V1.1.0 (2017-11), section 5.1
%   Data: 2017.12.25
%   Author: Zanyu, Chen   
%
%==========================================================================
    
function CRC_result = CRC_Check_NR_V2(data_CRC, CRC_case )
      
    CRC_tmp = Compute_CRC32_C_v2(data_CRC, CRC_case);
    if (CRC_tmp == 0)
        CRC_result = 1;
    else
        CRC_result = 0;
    end
end

