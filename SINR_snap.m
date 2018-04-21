%C:\Program Files\MATLAB\R2013b\cvx\MVDR


%%目前绘图结果不理想,有待研究原因
clear all;
clc

    snap=2:20:200;

monte_carloN=1;


sinr=zeros(1,length(snap));
sinr1=sinr;
sinr2=sinr;
sinr3=sinr;
sinr4=sinr;
ti=1;
     for N=snap %%fig=5
        snr_signal=10;
    for monte_carlo=1:monte_carloN
        q=4;
        sk=q-1;
        aa=[-60, -30, 30,0]*pi/180;
        M=7;
        m=M;
        n=[0:M-1]';
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
        snr=[30 30 30 snr_signal];
        [S0,Noise]=SNR(S,num,M,N,snr);
        A=zeros(M,q);
        for k=1:length(aa)
            A(:,k)=varray_line(aa(k),M,f(k),d);
        end
        X=A*S0+Noise;
        a=varray_line(0,M,1e9,d);% a=varray_line(doa,M,f,d)
%         a0=A(:,length(aa));% a=varray_line(doa,M,f,d)
        R=X*X'/N;
        Xjn=A(:,1:length(aa)-1)*S0(1:length(aa)-1,:)+Noise;
        Xs=A(:,length(aa))*S0(length(aa),:);
        Rs=Xs*Xs'/N;
        Rjn=Xjn*Xjn'/N;
        woo=inv(R)*a/(a'*inv(R)*a);
        %% 基于对角加载的MVDR，加载量10倍的噪声功率
        RDL=R+10*eye(M);
invRDL=inv(RDL);
w2=invRDL*a/(a'*invRDL*a);

%%%%%%%%%%%%%%%%%%%%%%%%%%method 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
as=(us*us'*SteerVec_sitai')';%a向信号子空间投影
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %   as= SteerVec_sitai;             
             
        

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
C=1;
dirta=3;
P=120;
g=1;
Kmax=999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
       
        ybsn=0.0001;
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
    w(:,K)=zeros(14,1);
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
         % jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)-hs(jw(:,K),Rxp,u(K,:),bi,ybsn,C)
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

w1=wVec';










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% proposed
pa=eye(M)-a*a'/(a'*a)+0.002*eye(M);
R3=pa*R*pa';
R3=(R3+R3')/2;
%% 降序取特征值特征向量
%需已知信号个数
[V,D] = eig(R3);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
ed=D_sort(1:q-1).*eye(q-1);%大特征值
us=V_sort(:,1:q-1);%大特征值对应的特征向量

B=inv(pa)*us;
%w2=-0.5*((-2*a'*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a-2)*inv(a'*a)*a+2*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a);

        


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sinr(ti)=sinr(ti)+10*log10((woo'*Rs*woo)/(woo'*Rjn*woo));
    sinr1(ti)=sinr1(ti)+10*log10((w1'*Rs*w1)/(w1'*Rjn*w1));
    sinr2(ti)=sinr2(ti)+10*log10((w2'*Rs*w2)/(w2'*Rjn*w2));
%     sinr3(ti)=sinr3(ti)+10*log10((w3'*Rs*w3)/(w3'*Rjn*w3));
%     sinr4(ti)=sinr4(ti)+10*log10((w4'*Rs*w4)/(w4'*Rjn*w4));
    end;
    sinr(ti)=sinr(ti)/monte_carloN;
    sinr1(ti)=sinr1(ti)/monte_carloN;
    sinr2(ti)=sinr2(ti)/monte_carloN;
%     sinr3(ti)=sinr3(ti)/monte_carloN;
%     sinr4(ti)=sinr4(ti)/monte_carloN;
    ti=ti+1;
end;

   
plot(snap,sinr,'k','LineWidth',1.0,'Marker','.');
hold on
plot(snap,sinr1,'b','LineWidth',1.0,'Marker','s');
hold on
grid on;
plot(snap,sinr2,'r','LineWidth',1.0,'Marker','p');
% plot(snap,sinr3,'g','LineWidth',1.0,'Marker','x');
% plot(snap,sinr4,'m','LineWidth',1.0,'Marker','o');
% grid on; hold on;
legend('MVDR','svm-mvdr','DL');
xlabel('Number of Snapshots');
ylabel('SINR/dB');
%% 
% axes('position', [0.17 0.7 0.29 0.18])
% plot(snap,sinr,'k','LineWidth',1.0,'Marker','.');
% hold on
% plot(snap,sinr3,'g','LineWidth',1.0,'Marker','x');
% grid on; hold on;
% plot(snap,sinr4,'m','LineWidth',1.0,'Marker','o');
% axis([100 300 36 38])
%%


