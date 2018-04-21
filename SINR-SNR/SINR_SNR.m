
% snrinput=-20:5:50;
% monte_carloN=100;
snrinput=-20:10:50;
monte_carloN=1;
sinr=zeros(1,length(snrinput));
sinr1=sinr;
sinr2=sinr;
ti=1;
for snr_signal=snrinput %%fig=4
        N=128;
    for monte_carlo=1:monte_carloN
        q=4;
        aa=[-60, 0, 30,1]*pi/180;
        M=5;
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
        w=inv(R)*a/(a'*inv(R)*a);
         %% 基于对角加载的MVDR，加载量10倍的噪声功率
RDL=R+10*eye(M);
invRDL=inv(RDL);
w1=invRDL*a/(a'*invRDL*a);
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
% B=A(:,2:sk);
B=inv(pa)*us;
w2=-0.5*((-2*a'*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a-2)*inv(a'*a)*a+2*B*inv(a'*a*B'*B-B'*a*a'*B)*B'*a);       
        
%%%%%%%%%%%%%%%%%%        
    sinr(ti)=sinr(ti)+10*log10((w'*Rs*w)/(w'*Rjn*w));
    sinr1(ti)=sinr1(ti)+10*log10((w1'*Rs*w1)/(w1'*Rjn*w1));
    sinr2(ti)=sinr2(ti)+10*log10((w2'*Rs*w2)/(w2'*Rjn*w2));
    end;
    sinr(ti)=sinr(ti)/monte_carloN;
    sinr1(ti)=sinr1(ti)/monte_carloN;
    sinr2(ti)=sinr2(ti)/monte_carloN;
    ti=ti+1;
end;

plot(snrinput,sinr,'k','LineWidth',1.0,'Marker','.');
hold on
plot(snrinput,sinr1,'b','LineWidth',1.0,'Marker','s');
hold on
grid on;
plot(snrinput,sinr2,'r','LineWidth',1.0,'Marker','p');

legend('MVDR','DL','PROPOSED');
xlabel('SNR/dB');
ylabel('SINR/dB');
