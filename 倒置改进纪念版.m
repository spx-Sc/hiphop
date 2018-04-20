% MVDR 波束形成的MATLAB仿真程序
% 包含经典MVDR
% Yujie Gu的协方差矩阵重构  Robust Adaptive Beamforming Based on Interference Covariance Matrix Reconstruction and Steering Vector Estimation
%208easy loop warning
% 自己提的一个小想法
clc; 
clear all; 
theta=[-2,60,20,-60];
ganrao=theta(2:size(theta,2));
num_j=size(ganrao,2);

snr=[-20,30,30,30];
m=10;%天线个数
a=array_line(0*pi/180,m);   %信号来向在这里规定？。。。。
sk=3;%信号个数
A=array_line(theta(1:sk)*pi/180,m);
N=256;
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
[S0,Noise]=SNR(S(1:sk,:),sk,m,N,snr(1:sk));   %产生真实传输的信号3*256 和各天线不同时刻产生的热噪音10*256
% %% 幅相误差
% bb=zeros(m-1,sk); ta=0.9; tb=1.1; tj=2;%相位
% tc=ta+(tb-ta)*rand(m-1,sk);%幅度
% bb=exp(j*tj*rand(m-1,sk)/180*pi);%相位
%   tb=tc.*bb;
%   Bg=[ones(1,sk); tb];
%   A=(A.*Bg);
% %%存在弗幅向误差
% bb=zeros(m-1,sk);  
% ta=0.9; tb=1.1; %幅度变化范围
% ta1=-2;tb1=2;%相位变化范围
% tc=ta+(tb-ta)*rand(m-1,sk);%幅度
% tc1=ta1+(tb1-ta1)*rand(m-1,sk);%相位
% bb=exp(1j*tc1/180*pi);%相位
%   tb=tc.*bb;
%   Bg=[ones(1,sk); tb];
%   A=(A.*Bg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=A*S0+Noise;          %接收到的信号，包含噪声、干扰、有用信息
% XI=A(:,2:sk)*S0(2:sk,:)+Noise; 
% RIN=XI*XI'/N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%常规MVDR
R=X*X'/N;
invR=inv(R);
w1=invR*a/(a'*invR*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基于对角加载的MVDR，加载量10倍的噪声功率
RDL=R+5*eye(m);  %对焦加载，因数就是5
invRDL=inv(RDL);
w2=invRDL*a/(a'*invRDL*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 特征子空间 ESB/应用了子空间投影理论
[V,D] = eig(R);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
us=V_sort(:,1:sk);
ed=D_sort(1:sk).*eye(sk);
a3=us*us'*a;%a向信号子空间投影
w3=invR*a3/(a3'*invR*a3);             %% beamforming
clear V D D_sort index V_sort us ed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% proposed
pa=eye(m)-a*a'/(a'*a)+0.002*eye(m);
R3=pa*R*pa';
R3=(R3+R3')/2;
%% 降序取特征值特征向量
%需已知信号个数
[V,D] = eig(R3);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
ed=D_sort(1:sk-1).*eye(sk-1);%大特征值
us=V_sort(:,1:sk-1);%大特征值对应的特征向量
Aq=inv(pa)*us;
%% CVX工具箱求解
cvx_begin
variable w4(m,1) complex;
minimize((w4'*w4))
% minimize(norm(w4))
 subject to
  w4'*a==1;
  w4'*inv(pa)*us==zeros(1,sk-1);
cvx_end


%%第五种解法
Kmax=9999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
        g=1; % 我们期望的反应
        ybsn=0.0001;
C=1;
dirta=2;
P=120;
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
        
        
      
%         for ff=1:60
%  %           SteerVec_sitai(ff,:) = SteerVecGen(TypeArray,sitai(ff),45*pi/180,M);
%                 
%         end
        as= SteerVec_sitai;
        
        
        
%         for i=1:60
%             aba(2*i-1,1)=real(as(i,:));  %aba is the a-ba; aba(120,2)
%             aba(2*i,1)=imag(as(i,:));
%             aba(2*i+119,1)=imag(as(i,:));
%             aba(2*i+120,1)=-real(as(i,:));
%         end
        for i=1:60
            denglou(:,i)=[real(as(i,:)') ; imag(as(i,:)')];
            denglou(:,i+60)=[imag(as(i,:)') ; -real(as(i,:)')];
        end
% 　　i=1;
%         aba1=[real(as(i,:)') ; imag(as(i,:)')];
%         aba2=[imag(as(i,:)') ; -real(as(i,:)')];
%         for i=2:60^^
%             aba1=[aba1 ; real(as(i,:)') ; imag(as(i,:)')];
%             aba2=[aba2 ; imag(as(i,:)') ; -real(as(i,:)')];
%         end
%         denglou=[aba1;aba2];
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
    while (  jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)-10000>=jhs(w(:,K),Rxp,u(K,:),bi,ybsn,C)   )
        yita(K)=yita(K)-0.1;                 % how fast the number is decreadsed in not sure
        w(:,K+1)=w(:,K)+yita(K)*(ws-w(:,K));
        for i=1:120    % warning!
            aba=denglou(:,i);
            u(K+1,i)=abs(di(i)-w(:,K+1)'*aba);
        end
          jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)
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


%% 闭式解
% B=A(:,2:sk);
% Aq=inv(pa)*us;
% w3=-0.5*((-2*a'*Aq*inv(a'*a*Aq'*Aq-Aq'*a*a'*Aq)*Aq'*a-2)*inv(a'*a)*a+2*Aq*inv(a'*a*Aq'*Aq-Aq'*a*a'*Aq)*Aq'*a);
% B=inv(m*Aq'*Aq-Aq'*a*a'*Aq);
% w3=(a/m)*(a'*Aq*B*Aq'*a+1)-Aq*B*Aq'*a;

%%
sita=[-90:90];                    %% 扫描方向范围
v=array_line(sita*pi/180,m);            %% 扫描方向矢量
B1=abs(w1'*v).^2; %MVER
B2=abs(w2'*v).^2; %
B3=abs(w3'*v).^2; %ESB
B4=abs(w4'*v).^2; %PROPOSED
B5=abs(w5'*v).^2; 
%%
figure;
plot(sita,10*log10(B1/max(B1)),'k:','LineWidth',1.0);
hold on
plot(sita,10*log10(B2/max(B2)),'r-','LineWidth',1.0);
plot(sita,10*log10(B3/max(B3)),'g-.','LineWidth',1.0);
plot(sita,10*log10(B4/max(B4)),'b--','LineWidth',1.0);
 plot(sita,10*log10(B5/max(B5)),'LineWidth',1.0);
B5=flipud(B5')';
%plot(,10*log10(B5/max(B5)),'black--','LineWidth',1.0);

 %SitaVec=linspace(-0.5*pi,0.5*pi,181).'
 
 
 % plot(SitaVec*180/pi,10*log10(B5/max(B5)),'b.-');
 
 
 plot(sita,10*log10(B5/max(B5)),'black*','LineWidth',1.0);
 
% hold on
% for i=1:sk  %画方位线
% plot([theta(i) theta(i)],[-90 0],'b:','LineWidth',1.0);
% end
xlabel('DOA/degree');ylabel('Beampattern/dB');
% legend('MVDR','DL','ESB','Proposed Method','extra','extraplus');
grid on
hold on
% axis([-90 90 -90 0]);
        for i=1:num_j
            yer=linspace(0,-70,201);
            xer=ones(201,1)*ganrao(i);
            plot(xer,yer);
        end


% pow=sum(abs(xs(1,:)).^2)/length(xs(1,:))
% 10*log10(pow)