function bit_2_img(A,H,W,CH)
%     filename="C:\Users\user\Desktop\PultoSDR\PLUTO\Linux_Share\0606_BIN_DATA.txt";
%     A=textread(filename,'%n','delimiter', ',');
    

    % A=textscan(filename,'%n','delimiter', ',');
    
    
    % filename2="0527_reDATA _1.txt";
    % K=textread(filename2,'%n','delimiter', ',');
    
    k=length(A)/8;
%     H=640;%要自己設定
%     W=480;%要自己設定
%     CH=3;%CH=3 彩色照片  CH=1 黑白照片
    i=1;
    B=[];
    B = zeros(1,k); %我加的
    while i<=k
        B(i)=strcat(string(A((i-1)*8+1)),string(A((i-1)*8+2)),string(A((i-1)*8+3)),string(A((i-1)*8+4)),string(A((i-1)*8+5)),string(A((i-1)*8+6)),string(A((i-1)*8+7)),string(A((i-1)*8+8)));
        i=i+1;
    end
    i=1;
    c=[];
    c = zeros(1,k); %我加的
    while i<=k
        c(i)=bin2dec(string(B(i)));
        i=i+1;
    end
    if CH==3
        i=1;
        c=[];
        c = zeros(1,k); %我加的
        while i<=k
            c(i)=bin2dec(string(B(i)));
            i=i+1;
        end
        D1=c(1:3:end);
        D2=c(2:3:end);
        D3=c(3:3:end);
    
        E1=reshape(D1,H,W)';
        E2=reshape(D2,H,W)';
        E3=reshape(D3,H,W)';
        F(:,:,1)=E3;
        F(:,:,2)=E2;
        F(:,:,3)=E1;
    elseif CH==1
        i=1;
        c=[];
        c = zeros(1,k); %我加的
        while i<=k
          c(i)=bin2dec(string(B(i)));
          i=i+1;
        end
        D1=c;
        E1=reshape(D1,W,H)';
        F=E1;
    
    end
    
    imshow(uint8(F));
end