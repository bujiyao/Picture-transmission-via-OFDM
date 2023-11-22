% =========================================================================
%   LTE modulation function
%   Support BPSK, QPSK, 16QAM, 64QAM.
%   References:
%   [1] 3GPP TS 36.211 V8.5.0 (2008-12), section 7.1.
%
%   Date: 2009.01.15
%   Author:  Fang, Yu-Chuan
%==========================================================================

function xk=modulation_LTE(bk,N_cpc)
type1=N_cpc;

switch type1
    case 1  %--- BPSK ---%
        tmp=-bk*2+1;
        xk = (tmp+i*tmp);  
        xk=xk/sqrt(2);
    case 2  %--- QPSK ---%
        xk=zeros(1,length(bk)/N_cpc);
        cnt_bk=1;
        for n=1:length(xk)
            xk(n)=(-bk(cnt_bk)*2+1);
            cnt_bk=cnt_bk+1;
            xk(n)=xk(n)+i*(-bk(cnt_bk)*2+1);
            cnt_bk=cnt_bk+1;
        end
        xk=xk/sqrt(2);
        
    case 4  %--- 16-QAM ---%
        xk=zeros(1,length(bk)/N_cpc);
        cnt_bk=1;
        for n=1:length(xk)
            %--- I channel -----%
            if bk(cnt_bk)==1 && bk(cnt_bk+2)==1
                xk(n)=xk(n)-3;
            elseif bk(cnt_bk)==1 && bk(cnt_bk+2)==0
                xk(n)=xk(n)-1;
            elseif bk(cnt_bk)==0 && bk(cnt_bk+2)==0
                xk(n)=xk(n)+1;
            else
                xk(n)=xk(n)+3;
            end
            %--- Q channel -----%
            if bk(cnt_bk+1)==1 && bk(cnt_bk+3)==1
                xk(n)=xk(n)-i*3;
            elseif bk(cnt_bk+1)==1 && bk(cnt_bk+3)==0
                xk(n)=xk(n)-i*1;
            elseif bk(cnt_bk+1)==0 && bk(cnt_bk+3)==0
                xk(n)=xk(n)+i*1;
            else
                xk(n)=xk(n)+i*3;
            end

            cnt_bk=cnt_bk+N_cpc;
        end
        
        xk=xk/sqrt(10);
        
    case 6  %--- 64-QAM ---%
        xk=zeros(1,length(bk)/N_cpc);
        cnt_bk=1;
        for n=1:length(xk)
            %--- I channel -----%
            if bk(cnt_bk)==1 && bk(cnt_bk+2)==1 && bk(cnt_bk+4)==1
                xk(n)=xk(n)-7;
            elseif bk(cnt_bk)==1 && bk(cnt_bk+2)==1 && bk(cnt_bk+4)==0
                xk(n)=xk(n)-5;
            elseif bk(cnt_bk)==1 && bk(cnt_bk+2)==0 & bk(cnt_bk+4)==0
                xk(n)=xk(n)-3;
            elseif bk(cnt_bk)==1 && bk(cnt_bk+2)==0 && bk(cnt_bk+4)==1
                xk(n)=xk(n)-1;
            elseif bk(cnt_bk)==0 && bk(cnt_bk+2)==0 && bk(cnt_bk+4)==1
                xk(n)=xk(n)+1;
            elseif bk(cnt_bk)==0 && bk(cnt_bk+2)==0 && bk(cnt_bk+4)==0
                xk(n)=xk(n)+3;
            elseif bk(cnt_bk)==0 && bk(cnt_bk+2)==1 && bk(cnt_bk+4)==0
                xk(n)=xk(n)+5;
            else
                xk(n)=xk(n)+7;
            end
            %--- Q channel -----%
            if bk(cnt_bk+1)==1 && bk(cnt_bk+3)==1 && bk(cnt_bk+5)==1
                xk(n)=xk(n)-i*7;
            elseif bk(cnt_bk+1)==1 && bk(cnt_bk+3)==1 && bk(cnt_bk+5)==0
                xk(n)=xk(n)-i*5;
            elseif bk(cnt_bk+1)==1 && bk(cnt_bk+3)==0 && bk(cnt_bk+5)==0
                xk(n)=xk(n)-i*3;
            elseif bk(cnt_bk+1)==1 && bk(cnt_bk+3)==0 && bk(cnt_bk+5)==1
                xk(n)=xk(n)-i*1;
            elseif bk(cnt_bk+1)==0 && bk(cnt_bk+3)==0 && bk(cnt_bk+5)==1
                xk(n)=xk(n)+i*1;
            elseif bk(cnt_bk+1)==0 && bk(cnt_bk+3)==0 && bk(cnt_bk+5)==0
               xk(n)=xk(n)+i*3;
            elseif bk(cnt_bk+1)==0 && bk(cnt_bk+3)==1 && bk(cnt_bk+5)==0
                xk(n)=xk(n)+i*5;
            else
                xk(n)=xk(n)+i*7;
            end
            
            cnt_bk=cnt_bk+N_cpc;
        end

        xk=xk/sqrt(42);
case 8   %--- 256-QAM ---%
        xk=zeros(1,length(bk)/N_cpc);
        cnt_bk=1;
        for n=1:length(xk)
            bit_0=(1-2*bk(cnt_bk));
            bit_2=(1-2*bk(cnt_bk+2));
            bit_4=(1-2*bk(cnt_bk+4));
            bit_6=(1-2*bk(cnt_bk+6));
            xk_real=bit_0*(8-bit_2*(4-bit_4*(2-bit_6)));

            bit_1=(1-2*bk(cnt_bk+1));
            bit_3=(1-2*bk(cnt_bk+3));
            bit_5=(1-2*bk(cnt_bk+5));
            bit_7=(1-2*bk(cnt_bk+7));
            xk_imag=j*bit_1*(8-bit_3*(4-bit_5*(2-bit_7)));
            cnt_bk=cnt_bk+8;

            xk(n)=xk_real+xk_imag;
        end
      
        xk=xk/sqrt(170);
        
    case 10   %--- 1024-QAM ---%
        xk=zeros(1,length(bk)/N_cpc);
        cnt_bk=1;
        for n=1:length(xk)
            bit_0=(1-2*bk(cnt_bk));
            bit_2=(1-2*bk(cnt_bk+2));
            bit_4=(1-2*bk(cnt_bk+4));
            bit_6=(1-2*bk(cnt_bk+6));
            bit_8=(1-2*bk(cnt_bk+8));
            xk_real=bit_0*(16-bit_2*(8-bit_4*(4-bit_6*(2-bit_8))));

            bit_1=(1-2*bk(cnt_bk+1));
            bit_3=(1-2*bk(cnt_bk+3));
            bit_5=(1-2*bk(cnt_bk+5));
            bit_7=(1-2*bk(cnt_bk+7));
            bit_9=(1-2*bk(cnt_bk+9));
            xk_imag=j*bit_1*(16-bit_3*(8-bit_5*(4-bit_7*(2-bit_9))));
            cnt_bk=cnt_bk+10;

            xk(n)=xk_real+xk_imag;
        end
        
        xk=xk/sqrt(682);

    otherwise
        disp('Unknown mapping')
end

