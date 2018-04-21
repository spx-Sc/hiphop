clear
clc
% snrinput=-20:5:50;
% monte_carloN=100;
%main 没有失配 并且角度做了修改,之前有个干扰来自0度,影响试验结果,改为25度

snrinput=-20:10:50;
monte_carloN=1;
sinr=zeros(1,length(snrinput));
sinr1=sinr;
sinr2=sinr;
ti=1;
%%%%%%%%%%%%%%%%%%%%%%w2-%%%%%%%%%%%%%%%%%%%%%%

Kmax=9999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
        g=1; % 我们期望的反应
        ybsn=0.0001;
C=1;
dirta=3;
P=120;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for snr_signal=snrinput %%fig=4
        N=128;
    for monte_carlo=1:monte_carloN
        q=4;
        sk=q-1;
        aa=[-60, -25, 20,0]*pi/180;
        M=5;
        m=M;
        d=3e8/1.2e9/2; 
        fc1=0.9;fc2=1;fc3=1.1;fc=1.2;
        f=[fc1,fc2,fc3,fc]*1e9;
        Ts=1/(4*max(f));
        
        fs=1/Ts;
        t=0:Ts:(N-1)*Ts;
        s1=exp(j*2*pi*6e6*t).*exp(j*2*pi*f(1)*t);
        s2=exp(j*2*pi*6e6*t).*exp(j*2*pi*f(2)*t);
        s3=exp(j*2*pi*6e4*t).*exp(j*2*pi*f(3)*t);
        s4=exp(-j*2*pi*6e6*t).*exp(j*2*pi*f(4)*t);
        S=[s1;s2;s3;s4];
        [k1,k2]=size(S);
        num=k1;
        snr=[25 25 25 snr_signal];
        [S0,Noise]=SNR(S,num,M,N,snr);
        A=zeros(M,q);
        for k=1:length(aa)
            A(:,k)=varray_line(aa(k),M,f(k),d);
        end
        X=A*S0+Noise;
        a=varray_line(0,M,1e9,d);% a=varray_line(doa,M,f,d)
        R=X*X'/N;
        Xjn=A(:,1:length(aa)-1)*S0(1:length(aa)-1,:)+Noise;
        Xs=A(:,length(aa))*S0(length(aa),:);
        Rs=Xs*Xs'/N;
        Rjn=Xjn*Xjn'/N;
        w0=inv(R)*a/(a'*inv(R)*a);
         %% 基于对角加载的MVDR，加载量10倍的噪声功率
RDL=R+10*eye(M);
invRDL=inv(RDL);
w1=invRDL*a/(a'*invRDL*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% proposed
% pa=eye(M)-a*a'/(a'*a)+0.002*eye(M);
% R3=pa*R*pa';
% R3=(R3+R3')/2;
% %% 降序取特征值特征向量
% %需已知信号个数
% [V,D] = eig(R3);
% [D_sort,index] = sort(diag(D),'descend');
% V_sort = V(:,index);
% ed=D_sort(1:q-1).*eye(q-1);%大特征值
% us=V_sort(:,1:q-1);%大特征值对应的特征向量
% %% 用幂法求解信号子空间，需已知信号个数
% for i=1:sk-1  
% [ed(i,1), u(:,i)] = powermethod(R3);%用幂法求解信号子空间，需已知信号个数
% u(:,i)=u(:,i)/norm(u(:,i));
% R3=R3-ed(i,1)*u(:,i)*u(:,i)';
% end
% us=u;
% %% CVX工具箱求解
% cvx_begin
% variable w4(m,1) complex;
% minimize((w4'*w4))
%  subject to
%   w4'*a==1;
%   w4'*inv(pa)*us==zeros(1,sk-1);
% cvx_end
%% 闭式解
% % B=A(:,2:sk);
% B=inv(pa)*us;
% w2=-0.5*((-2*a'*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a-2)*inv(a'*a)*a+2*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a);       
%         
%%%%%%%%%%%new w2 our method

%%第五种解法
        sitai=linspace(-90,90,60);
        SteerVec_sitai=array_linerr(sitai,m)';
        
             for ff=1:60
             SteerVec_sitai(ff,:) = SteerVecGen('ULA',sitai(ff)*pi/180,45*pi/180,m);
 %               SteerVec_sitai(ff,:)=array_linerr(sitai(ff)*pi/180,m)';
             end   
        
        
      
            %%%%%%%%%%%%%%%%%%%%%%esb分解实验%%%%%%%%%%%%%%%%
[V,D] = eig(R);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
us=V_sort(:,1:sk);
ed=D_sort(1:sk).*eye(sk);
as=(us*us'*SteerVec_sitai')';%a向信号子空间投影
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

        for i=1:60
            denglou(:,i)=[real(as(i,:)') ; imag(as(i,:)')];
            denglou(:,i+60)=[imag(as(i,:)') ; -real(as(i,:)')];
        end
Rx=R;      

Rxp=[real(Rx),-imag(Rx);imag(Rx),real(Rx)];   
         
  
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
      
u=(di);   %u's tructure : u(k,i) ,it should be in the outer of loop
for K=1:Kmax
    w(:,K)=zeros(2*m,1);
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
    ws=(Rxp+denglou*Df*denglou')\(denglou*Df*di');
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
    
    if ( jhs(w(:,kjia),Rxp,u(kjia,:),bi,ybsn,C)- jhs(w(:,kjia-1),Rxp,u(kjia-1,:),bi,ybsn,C))<minchaju;
        break;
    end
end
for i=1:m
    wVec(i)=w(i,K+1)+w(i+m,K+1)*j;
    wjia(i)=w(i,K+1)^2+(w(i+m,K+1)*j)^2;
end
%wVec=wVec';

w2=wVec;

if size(w2,1)==1
    w2=w2';
end

%%%%%%%%%%%%%%%%%%        
    sinr(ti)=sinr(ti)+10*log10((w0'*Rs*w0)/(w0'*Rjn*w0));
    sinr1(ti)=sinr1(ti)+10*log10((w1'*Rs*w1)/(w1'*Rjn*w1));
    sinr2(ti)=sinr2(ti)+10*log10((w2'*Rs*w2)/(w2'*Rjn*w2));
    end;
    sinr(ti)=sinr(ti)/monte_carloN;
    sinr1(ti)=sinr1(ti)/monte_carloN;
    sinr2(ti)=sinr2(ti)/monte_carloN;
    ti=ti+1;
end;
figure;
plot(snrinput,sinr,'k','LineWidth',1.0,'Marker','.');
hold on
plot(snrinput,sinr1,'b','LineWidth',1.0,'Marker','s');
hold on
grid on;
plot(snrinput,sinr2,'r','LineWidth',1.0,'Marker','p');
hold off;
legend('MVDR','DL','SVM-ESB');
xlabel('SNR/dB');
ylabel('SINR/dB');
