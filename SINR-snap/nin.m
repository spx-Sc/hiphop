function a=nin()
%����esb��ʱ����.m
a=1;
theta=[3,-30,20];
ganrao=theta(2:size(theta,2));
num_j=size(ganrao,2);
snr=[0,30,30,30];
m=10;%���߸���
a=array_line(0*pi/180,m);   %�ź�����������涨����������
sk=3;%�źŸ���
A=array_line(theta(1:sk)*pi/180,m);
N=256;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kmax=9999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
        g=1; % ���������ķ�Ӧ
        ybsn=0.0001;
C=1;
dirta=5;
P=120;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







fc1=0.9194;
fc2=1;
fc3=0.8387;%520*10^6;
Ts=1/(4*fc2);
t=0:Ts:(N-1)*Ts;
s1=1*sin(2*pi*fc1*t);
s1=hilbert(s1);
s2=sin(2*pi*fc2*t+5*cos(80*t));
s2=hilbert(s2);
s3=1*sin(2*pi*fc3*t+5*cos(10*t));
s3=hilbert(s3);
% fc=1;
% N_sample=160;
% N1=100;
% B=1;
% [st,ss]=mbpsk(fc,N_sample,N1,B);
% s0=ss(1:N);
signal = generatesignal(1,1,10);    
s0=signal(1:N);
S=[s0;s1;s2;s3];
[S0,Noise]=SNR(S(1:sk,:),sk,m,N,snr(1:sk));   %������ʵ������ź�3*256 �͸����߲�ͬʱ�̲�����������10*256
% %% �������
% bb=zeros(m-1,sk); ta=0.9; tb=1.1; tj=2;%��λ
% tc=ta+(tb-ta)*rand(m-1,sk);%����
% bb=exp(j*tj*rand(m-1,sk)/180*pi);%��λ
%   tb=tc.*bb;
%   Bg=[ones(1,sk); tb];
%   A=(A.*Bg);
% %%���ڸ��������
% bb=zeros(m-1,sk);  
% ta=0.9; tb=1.1; %���ȱ仯��Χ
% ta1=-2;tb1=2;%��λ�仯��Χ
% tc=ta+(tb-ta)*rand(m-1,sk);%����
% tc1=ta1+(tb1-ta1)*rand(m-1,sk);%��λ
% bb=exp(1j*tc1/180*pi);%��λ
%   tb=tc.*bb;
%   Bg=[ones(1,sk); tb];
%   A=(A.*Bg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=A*S0+Noise;          %���յ����źţ��������������š�������Ϣ
% XI=A(:,2:sk)*S0(2:sk,:)+Noise; 
% RIN=XI*XI'/N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����MVDR
R=X*X'/N;
invR=inv(R);
w1=invR*a/(a'*invR*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ڶԽǼ��ص�MVDR��������10������������
RDL=R+5*eye(m);  %�Խ����أ���������5
invRDL=inv(RDL);
w2=invRDL*a/(a'*invRDL*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����ӿռ� ESB/Ӧ�����ӿռ�ͶӰ����
[V,D] = eig(R);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
us=V_sort(:,1:sk);
ed=D_sort(1:sk).*eye(sk);
a3=us*us'*a;%a���ź��ӿռ�ͶӰ
w3=invR*a3/(a3'*invR*a3);             %% beamforming
clear V D D_sort index V_sort us ed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% proposed

%% ����ȡ����ֵ��������
%����֪�źŸ���


%%�����ֽⷨ
%         sitai(1:40)=linspace(-90,90,40);
%         sitai(41:80)=linspace(-90,-5,40);
%         sitai(81:120)=linspace(5,90,40);
%         for ff=1:120
%         sitai(1:20)=linspace(-90,90,20);
%         sitai(21:40)=linspace(-90,-5,20);
%         sitai(41:60)=linspace(5,90,20);  
        sitai=linspace(-90,90,60);
        SteerVec_sitai=array_linerr(sitai,m)';
        
             for ff=1:60
             SteerVec_sitai(ff,:) = SteerVecGen('ULA',sitai(ff)*pi/180,45*pi/180,m);
 %               SteerVec_sitai(ff,:)=array_linerr(sitai(ff)*pi/180,m)';
             end   
        
        
      
            %%%%%%%%%%%%%%%%%%%%%%esb�ֽ�ʵ��%%%%%%%%%%%%%%%%
[V,D] = eig(R);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
us=V_sort(:,1:sk);
ed=D_sort(1:sk).*eye(sk);
%as=(us*us'*SteerVec_sitai')';%a���ź��ӿռ�ͶӰ
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        as= SteerVec_sitai;             
             
        

        for i=1:60
            denglou(:,i)=[real(as(i,:)') ; imag(as(i,:)')];
            denglou(:,i+60)=[imag(as(i,:)') ; -real(as(i,:)')];
        end
Rx=R;      

Rxp=[real(Rx),-imag(Rx);imag(Rx),real(Rx)];   %Rxp!=Rx!!!!!!!!!!!!!
         
        
        
%         for i=1:120
%             if abs(sitai(i))>dirta
%                 di(i)=0;
%             else
%                 di=g;
%             end
%         end
        for i=1:60
            if abs(sitai(i))>dirta
                di(i)=0;
            else
                di(i)=g;
            end
            di(i)=real(di(i));
            if abs(sitai(i))>dirta
                di(i+60)=0;
            else
                di(i+60)=g;
            end
            di(i+60)=imag(di(i+60));
        end
       % Rx=RxMat;
        
        
%         u(1,:)=(di);   %u's tructure : u(k,i) ,it should be in the outer of loop
%         for K=1:Kmax
%             w(:,K)=0;
%             for i=1:120
%                 if u(K,i)<ybsn
%                     fi(i)=0;
%                 else
%                     fi(i)=C/u(K,i);
%                 end
%                 Df(i,i)=fi(i);
%             end
u=(di);   %u's tructure : u(k,i) ,it should be in the outer of loop
for K=1:Kmax
    w(:,K)=zeros(20,1);
    for i=1:120                       %using method m=1 another place used
        if u(K,i)<ybsn
            fi(i)=0;
        else
           % fi(i)=C/u(K,i);   j=1
           fi(i)=2*C*(u(K,i)-ybsn)/u(K,i);    %m=2version
        end
        Df(i,i)=fi(i);
    end
    %%%%%
    ws=(Rxp+denglou*Df*denglou')\(denglou*Df*di');%%%the problem still exist.continue!
    yita(K)=1;
    %%%
    
    w(:,K+1)=w(:,K)+yita(K)*(ws-w(:,K));
    for i=1:120    % warning!
        aba=denglou(:,i);
        u(K+1,i)=abs(di(i)-w(:,K+1)'*aba);
    end
    while (  jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)-0.5>=jhs(w(:,K),Rxp,u(K,:),bi,ybsn,C)   )
        yita(K)=yita(K)-0.05;                 % how fast the number is decreadsed in not sure
        w(:,K+1)=w(:,K)+yita(K)*(ws-w(:,K));
        for i=1:120    % warning!
            aba=denglou(:,i);
            u(K+1,i)=abs(di(i)-w(:,K+1)'*aba);
        end
          jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)-jhs(w(:,K),Rxp,u(K,:),bi,ybsn,C)
    end
    % form now, the value of k is virtually increased.
    kjia=K+1;
    for i=1:120    % warning!
        aba=denglou(:,i);
        u(kjia,i)=abs(di(i)-w(:,kjia)'*aba);
    end
    %%%%%%%%
    for i=1:120     %using method m=1   another place used
        if u(kjia,i)<ybsn
            fi(i)=0;
        else
            fi(i)=C/u(kjia,i);
        end
    end
    %%%%%%%%%
%     jhs(w(kjia+1),Rxp,u(kjia+1,:),bi,ybsn,C)
%     jhs(w(:,kjia),Rxp,u(kjia,:),bi,ybsn,C)
   % if ((jhs(w(kjia+1))-jhs(w(:,kjia)))<minchaju )
    if ( jhs(w(:,kjia),Rxp,u(kjia,:),bi,ybsn,C)- jhs(w(:,kjia-1),Rxp,u(kjia-1,:),bi,ybsn,C))<minchaju;
        break;
    end
end
for i=1:10
    wVec(i)=w(i,K+1)+w(i+10,K+1)*j;
    wjia(i)=w(i,K+1)^2+(w(i+10,K+1)*j)^2;
end
wVec=wVec';
w5=wVec;

sita=[-90:90];                    %% ɨ�跽��Χ
v=array_line(sita*pi/180,m);            %% ɨ�跽��ʸ��

B5=abs(w5'*v).^2; 
B5=flipud(B5')';
 plot(sita,10*log10(B5/max(B5)),'cyan-.','LineWidth',1.0);