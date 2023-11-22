function [CellID,PSS,SSS] = SS(N1,N2)

% N1=0;N2=0;
CellID=3*N1+N2;

n1=[0:335];
n2=[0:2];
% % % % PSS
% % % % %  NR (m-seq)
x=[0,1,1,0,1,1,1];  %%% x1~x7=[[0,1,1,0,1,1,1]]
c=7;
for i = 1:121          %%  x8~128
    x(i+c)=mod(x(i+4)+x(i),2);
end

% % % % %  NR (m-seq) index

n=[0:126];
for i=0:2
m(i+1,:)=mod((n+43*i),127);
end
m=m+1;  %% index start 

% % % % NR PSS
for i=1:3
pss(i,:)=1-2*x(m(i,:));
end

switch  N2
    case 0
    PSS=pss(1,:);
    case 1
    PSS=pss(2,:);
    case 2
    PSS=pss(3,:);
end
% % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % NR   SSS
% % % % % % % % % % % % % % % % % % % 

% % % % % % %   m1
m1=mod(n1,112);
% % % % % m0
for i= 0:2
m0(i+1,:)=15*floor(n1/112)+5*i;
end





% % % % % % % 
x0=[1,0,0,0,0,0,0];  %%% x1~x7=[[0 0 0 0 0 0  1]]
c=7;
for i = 1:121          %%  x8~128
    x0(i+c)=mod(x0(i+4)+x0(i),2);
end
% % % % % 
% % % % % % % 
x1=[1,0,0,0,0,0,0];  %%% x1~x7=[[0 0 0 0 0 0  1]]
c=7;
for i = 1:121          %%  x8~128
    x1(i+c)=mod(x1(i+1)+x1(i),2);
end
% % % % % 
% % % % NR SSS
for k =1:3
for i=1:336
sss(i,:,k)=(1-2*x0(mod(n+m0(k,i),127)+1)).*(1-2*x1(mod(n+m1(k),127)+1));
end
end
i=N1+1;k=N2+1;
SSS=sss(i,:,k);




end


