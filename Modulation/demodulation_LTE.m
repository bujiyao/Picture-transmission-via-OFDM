% =========================================================================
%   LTE demodulation function
%   Output is the signal point on the constellation by the hard decision.
%   Support BPSK, QPSK, 16QAM, 64QAM.
%   References:
%   [1] 3GPP TS 36.211 V8.5.0 (2008-12), section 7.1.
%
%   Date: 2008.01.15
%   Author: Fang, Yu-Chuan
%==========================================================================

function yk=demodulation_LTE(xk,N_cpc)

type1 = N_cpc;
yk = zeros(size(xk));
Num_xk = length(xk);
switch type1
    case 1  %--- BPSK ---%
        re_xk = real(xk);
        im_xk = imag(xk);
        for n = 1:Num_xk
            if im_xk(n) >= -re_xk(n)
                yk(n) = 1+i;
            else
                yk(n) = -1-i;
            end
        end
         %yk=yk/sqrt(2);
    case 2  %--- QPSK ---%
        re_xk = real(xk);
        im_xk = imag(xk);
        re_yk = (re_xk>=0)*2.0-1;
        im_yk = (im_xk>=0)*2.0-1;
        yk = re_yk + i*im_yk;
%         yk=yk/sqrt(2);
 case 4  %--- 16-QAM ---%
        xk=xk*sqrt(10);
        re_xk = real(xk);
        im_xk = imag(xk);
        for n = 1:Num_xk
            if re_xk(n) >= 2
                yk(n) = 3;
            elseif (re_xk(n) < 2 && (re_xk(n) >= 0))
                yk(n) = 1;
            elseif (re_xk(n) < 0 && (re_xk(n) >= (-2)))
                yk(n) = -1;
            else
                yk(n) = -3;
            end
            if im_xk(n) >= 2
                yk(n) = yk(n)+i*3;
            elseif (im_xk(n) < 2 && (im_xk(n) >= 0))
                yk(n) = yk(n)+i;
            elseif (im_xk(n) < 0 && (im_xk(n) >= (-2)))
                yk(n) = yk(n)-i;
            else
                yk(n) = yk(n)-i*3;
            end
        end
    case 6  %--- 64-QAM ---%
        re_xk=real(xk);
        im_xk=imag(xk);
        for n=1:Num_xk
            if re_xk(n) >= (1.25-(5/16))
                yk(n) = 7;
            elseif (re_xk(n) < (1.25-(5/16))) && (re_xk(n) >=(1.25-(5/16)*2))
                yk(n) = 5;
            elseif (re_xk(n) < (1.25-(5/16)*2)) && (re_xk(n) >= (1.25-(5/16)*3))
                yk(n) = 3;
            elseif (re_xk(n) < (1.25-(5/16)*3)) && (re_xk(n) >= (1.25-(5/16)*4))
                yk(n) = 1;
            elseif (re_xk(n) < (1.25-(5/16)*4)) && (re_xk(n) >= (1.25-(5/16)*5))
                yk(n) = -1;
            elseif (re_xk(n) < (1.25-(5/16)*5)) && (re_xk(n) >= (1.25-(5/16)*6))
                yk(n) = -3;
            elseif (re_xk(n) < (1.25-(5/16)*6)) && (re_xk(n) >= (1.25-(5/16)*7))
                yk(n) = -5;
            else
                yk(n) = -7;
            end

            if im_xk(n) >= (1.25-(5/16))
                yk(n) = yk(n)+i*7;
            elseif (im_xk(n) < (1.25-(5/16))) && (im_xk(n) >=(1.25-(5/16)*2))
                yk(n) = yk(n)+i*5;
            elseif (im_xk(n) < (1.25-(5/16)*2)) && (im_xk(n) >= (1.25-(5/16)*3))
                yk(n) = yk(n)+i*3;
            elseif (im_xk(n) < (1.25-(5/16)*3)) && (im_xk(n) >= (1.25-(5/16)*4))
                yk(n) = yk(n)+i*1;
            elseif (im_xk(n) < (1.25-(5/16)*4)) && (im_xk(n) >= (1.25-(5/16)*5))
                yk(n) = yk(n)-i;
            elseif (im_xk(n) < (1.25-(5/16)*5)) && (im_xk(n) >= (1.25-(5/16)*6))
                yk(n) = yk(n)-3*i;
            elseif (im_xk(n) < (1.25-(5/16)*6)) && (im_xk(n) >= (1.25-(5/16)*7))
                yk(n) = yk(n)-i*5;
            else
                yk(n) = yk(n)-i*7;
            end
        end
%         yk=yk/sqrt(42);


 case 8  %--- 256-QAM ---%
        xk=xk*sqrt(170);
        re_xk = real(xk);
        im_xk = imag(xk);
        for n = 1:Num_xk
            if re_xk(n) >= 14
                yk(n) = 15;
            elseif (re_xk(n) < 14 && (re_xk(n) >= 12))
                yk(n) = 13;
            elseif (re_xk(n) < 12 && (re_xk(n) >= 10))
                yk(n) = 11;
                elseif (re_xk(n) < 10 && (re_xk(n) >= 8))
                yk(n) = 9;
                elseif (re_xk(n) < 8 && (re_xk(n) >= 6))
                yk(n) = 7;
                elseif (re_xk(n) < 6 && (re_xk(n) >= 4))
                yk(n) = 5;
                elseif (re_xk(n) < 4 && (re_xk(n) >= 2))
                yk(n) = 3;              
                elseif (re_xk(n) < 2 && (re_xk(n) >= 0))
                yk(n) = 1;
                elseif (re_xk(n) < 0 && (re_xk(n) >= (-2)))
                yk(n) = (-1);      
                elseif (re_xk(n) < (-2) && (re_xk(n) >= (-4)))
                yk(n) = (-3); 
                elseif (re_xk(n) < (-4) && (re_xk(n) >= (-6)))
                yk(n) = (-5);
                 elseif (re_xk(n) < (-6) && (re_xk(n) >= (-8)))
                yk(n) = (-7);
                elseif (re_xk(n) < (-8) && (re_xk(n) >= (-10)))
                yk(n) = (-9);
                elseif (re_xk(n) < (-10) && (re_xk(n) >= (-12)))
                yk(n) = (-11);
                 elseif (re_xk(n) < (-12) && (re_xk(n) >= (-14)))
                yk(n) = (-13);
            else
                yk(n) = -15;
            end
            if im_xk(n) >= 14
                yk(n) = yk(n)+i*15;
            elseif (im_xk(n) < 14 && (im_xk(n) >= 12))
                yk(n) = yk(n)+i*13;
            elseif (im_xk(n) < 12 && (im_xk(n) >= 10))
                yk(n) = yk(n)+i*11;
            elseif (im_xk(n) < 10 && (im_xk(n) >= 8))
                yk(n) = yk(n)+i*9;
            elseif (im_xk(n) < 8 && (im_xk(n) >= 6))   
                yk(n) = yk(n)+i*7;
            elseif (im_xk(n) < 6 && (im_xk(n) >= 4))
                yk(n) = yk(n)+i*5;
            elseif (im_xk(n) < 4 && (im_xk(n) >= 2))
                yk(n) = yk(n)+i*3;
            elseif (im_xk(n) < 2 && (im_xk(n) >= 0))
                yk(n) = yk(n)+i*1;
            elseif (im_xk(n) < 0 && (im_xk(n) >= (-2)))
                yk(n) = yk(n)-i*(1);
            elseif (im_xk(n) < (-2) && (im_xk(n) >= (-4)))
                yk(n) = yk(n)-i*(3);
            elseif (im_xk(n) < (-4) && (im_xk(n) >= (-6)))
                yk(n) = yk(n)-i*(5);
            elseif (im_xk(n) < (-6) && (im_xk(n) >= (-8)))
                yk(n) = yk(n)-i*(7);
            elseif (im_xk(n) < (-8) && (im_xk(n) >= (-10)))
                yk(n) = yk(n)-i*(9);
             elseif (im_xk(n) < (-10) && (im_xk(n) >= (-12)))
                yk(n) = yk(n)-i*(11);
                elseif (im_xk(n) < (-12) && (im_xk(n) >= (-14)))
                yk(n) = yk(n)-i*(13);
            else
                yk(n) = yk(n)-i*(15);
            end
        end
        
        
        
   case 10  %--- 1024-QAM ---%
        xk=xk*sqrt(682);
        re_xk = real(xk);
        im_xk = imag(xk);
        for n = 1:Num_xk
            if re_xk(n) >= 30
                yk(n) = 31;
            elseif (re_xk(n) < 30 && (re_xk(n) >= 28))
                yk(n) = 29;
            elseif (re_xk(n) < 28 && (re_xk(n) >= 26))
                yk(n) = 27;
                elseif (re_xk(n) < 26 && (re_xk(n) >= 24))
                yk(n) = 25;
                elseif (re_xk(n) < 24 && (re_xk(n) >= 22))
                yk(n) = 23;
                elseif (re_xk(n) < 22 && (re_xk(n) >= 20))
                yk(n) = 21;
                elseif (re_xk(n) < 20 && (re_xk(n) >= 18))
                yk(n) = 19;              
                elseif (re_xk(n) < 18 && (re_xk(n) >= 16))
                yk(n) = 17;
                elseif (re_xk(n) < 16 && (re_xk(n) >= 14))
                yk(n) = 15;
                elseif (re_xk(n) < 14 && (re_xk(n) >= 12))
                yk(n) = 13;
                elseif (re_xk(n) < 12 && (re_xk(n) >= 10))
                yk(n) = 11;
                elseif (re_xk(n) < 10 && (re_xk(n) >= 8))
                yk(n) = 9;
                elseif (re_xk(n) < 8 && (re_xk(n) >= 6))
                yk(n) = 7;              
                elseif (re_xk(n) < 6 && (re_xk(n) >= 4))
                yk(n) = 5;
                elseif (re_xk(n) < 4 && (re_xk(n) >= 2))
                yk(n) = 3;
                elseif (re_xk(n) < 2 && (re_xk(n) >= 0))
                yk(n) = 1;              
                elseif (re_xk(n) < 0 && (re_xk(n) >= (-2)))
                yk(n) = (-1);      
                elseif (re_xk(n) < (-2) && (re_xk(n) >= (-4)))
                yk(n) = (-3); 
                elseif (re_xk(n) < (-4) && (re_xk(n) >= (-6)))
                yk(n) = (-5);
                 elseif (re_xk(n) < (-6) && (re_xk(n) >= (-8)))
                yk(n) = (-7);
                elseif (re_xk(n) < (-8) && (re_xk(n) >= (-10)))
                yk(n) = (-9);
                elseif (re_xk(n) < (-10) && (re_xk(n) >= (-12)))
                yk(n) = (-11);
                 elseif (re_xk(n) < (-12) && (re_xk(n) >= (-14)))
                yk(n) = (-13);
                elseif (re_xk(n) < (-14) && (re_xk(n) >= (-16)))
                yk(n) = (-15); 
                elseif (re_xk(n) < (-16) && (re_xk(n) >= (-18)))
                yk(n) = (-17);
                 elseif (re_xk(n) < (-18) && (re_xk(n) >= (-20)))
                yk(n) = (-19);
                elseif (re_xk(n) < (-20) && (re_xk(n) >= (-22)))
                yk(n) = (-21);
                elseif (re_xk(n) < (-22) && (re_xk(n) >= (-24)))
                yk(n) = (-23);
                 elseif (re_xk(n) < (-24) && (re_xk(n) >= (-26)))
                yk(n) = (-25);
                elseif (re_xk(n) < (-26) && (re_xk(n) >= (-28)))
                yk(n) = (-27); 
                elseif (re_xk(n) < (-28) && (re_xk(n) >= (-30)))
                yk(n) = (-29);
  
            else
                yk(n) = -31;
            end
            if im_xk(n) >= 30
                yk(n) = yk(n)+i*31;
            elseif (im_xk(n) < 30 && (im_xk(n) >= 28))
                yk(n) = yk(n)+i*29;
            elseif (im_xk(n) < 28 && (im_xk(n) >= 26))
                yk(n) = yk(n)+i*27;
            elseif (im_xk(n) < 26 && (im_xk(n) >= 24))
                yk(n) = yk(n)+i*25;
            elseif (im_xk(n) < 24 && (im_xk(n) >= 22))   
                yk(n) = yk(n)+i*23;
            elseif (im_xk(n) < 22 && (im_xk(n) >= 20))
                yk(n) = yk(n)+i*21;
            elseif (im_xk(n) < 20 && (im_xk(n) >= 18))
                yk(n) = yk(n)+i*19;
            elseif (im_xk(n) < 18 && (im_xk(n) >= 16))
                yk(n) = yk(n)+i*17;
            elseif (im_xk(n) < 16 && (im_xk(n) >= 14))
                yk(n) = yk(n)+i*15;
            elseif (im_xk(n) < 14 && (im_xk(n) >= 12))
                yk(n) = yk(n)+i*13;
            elseif (im_xk(n) < 12 && (im_xk(n) >= 10))
                yk(n) = yk(n)+i*11;
            elseif (im_xk(n) < 10 && (im_xk(n) >= 8))   
                yk(n) = yk(n)+i*9;
                elseif (im_xk(n) < 8 && (im_xk(n) >= 6))   
                yk(n) = yk(n)+i*7;
            elseif (im_xk(n) < 6 && (im_xk(n) >= 4))
                yk(n) = yk(n)+i*5;
            elseif (im_xk(n) < 4 && (im_xk(n) >= 2))
                yk(n) = yk(n)+i*3;
            elseif (im_xk(n) < 2 && (im_xk(n) >= 0))
                yk(n) = yk(n)+i*1;    
                
                
            elseif (im_xk(n) < 0 && (im_xk(n) >= (-2)))
                yk(n) = yk(n)-i*(1);
            elseif (im_xk(n) < (-2) && (im_xk(n) >= (-4)))
                yk(n) = yk(n)-i*(3);
            elseif (im_xk(n) < (-4) && (im_xk(n) >= (-6)))
                yk(n) = yk(n)-i*(5);
            elseif (im_xk(n) < (-6) && (im_xk(n) >= (-8)))
                yk(n) = yk(n)-i*(7);
            elseif (im_xk(n) < (-8) && (im_xk(n) >= (-10)))
                yk(n) = yk(n)-i*(9);
             elseif (im_xk(n) < (-10) && (im_xk(n) >= (-12)))
                yk(n) = yk(n)-i*(11);
            elseif (im_xk(n) < (-12) && (im_xk(n) >= (-14)))
                yk(n) = yk(n)-i*(13);
            elseif (im_xk(n) < (-14) && (im_xk(n) >= (-16)))
                yk(n) = yk(n)-i*(15);
            elseif (im_xk(n) < (-16) && (im_xk(n) >= (-18)))
                yk(n) = yk(n)-i*(17);
            elseif (im_xk(n) < (-18) && (im_xk(n) >= (-20)))
                yk(n) = yk(n)-i*(19);
            elseif (im_xk(n) < (-20) && (im_xk(n) >= (-22)))
                yk(n) = yk(n)-i*(21);
             elseif (im_xk(n) < (-22) && (im_xk(n) >= (-24)))
                yk(n) = yk(n)-i*(23);
            elseif (im_xk(n) < (-24) && (im_xk(n) >= (-26)))
                yk(n) = yk(n)-i*(25);
                elseif (im_xk(n) < (-26) && (im_xk(n) >= (-28)))
                yk(n) = yk(n)-i*(27);
                elseif (im_xk(n) < (-28) && (im_xk(n) >= (-30)))
                yk(n) = yk(n)-i*(29);
                
            else
                yk(n) = yk(n)-i*(31);
            end
        end      
%  yk=yk/sqrt(10);
    otherwise
        disp('Unknown mapping')
end

