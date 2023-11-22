function [r] = CSIRS(n,l,Num_Used_Carrier)
    N=14; %number of symbol per slot
    % n=  1 ;   %is the slot number within a radio frame
    % l=  5 ;   %is the OFDM symbol within a slot
    nID= 0;    %
    CSI_C_init=mod(((2^10)*(N*n+l+1)*(2*nID+1)+nID),2^31);
    
    x1=[1,zeros(1,30)];
    x=31;
    for n=1:10000
        x1(n+x)=mod(x1(n+3)+x1(n),2);
    end
    xx=[0];
    cont=1;
    c1=log2(CSI_C_init);
    for i=30:-1:0
        if (floor(log2(CSI_C_init))-i ) == 0
            xx(cont)=[i];
           CSI_C_init= CSI_C_init-2^xx(cont);
            cont=cont+1;
        end
    end
    xx=xx+1;
    x2=[zeros(1,31)];
    for i=1:numel(xx)
    x2(xx(i))=1;
    end
    x=31;
    for n=1:10000
        x2(n+x)=mod((x2(n+3)+x2(n+2)+x2(n+1)+x2(n)),2);
    end
    
    for n= 1:Num_Used_Carrier*2+1
        C(n)=mod(x1(n+1600)+x2(n+1600),2);
    end
    for m=1:Num_Used_Carrier
        r(m)=((1/sqrt(2))*(1-2*C(2*m)))+j*(1/sqrt(2))*(1-2*C(2*m+1));
    end
end