% MVDR �����γɵ�MATLAB�������
% ��������MVDR
% Yujie Gu��Э��������ع�  Robust Adaptive Beamforming Based on Interference Covariance Matrix Reconstruction and Steering Vector Estimation

% �Լ����һ��С�뷨
clc; 
% close all
clear all; 
theta=[2,-30,60,20];
snr=[-15,30,30,40];
sk=3;%�źŸ���
m=10;%���߸���
n=[0:m-1]'; %n
A=array_line(theta(1:sk)*pi/180,m);
N=256;
fc1=0.9194;
fc2=1;
fc3=0.8387;%520*10^6;
d=0.15;
Ts=1/(4*fc2);
t=0:Ts:(N-1)*Ts;
s1=1*sin(2*pi*fc1*t);
s1=hilbert(s1);
s2=sin(2*pi*fc2*t+5*cos(80*t));
s2=hilbert(s2);
s3=1*sin(2*pi*fc3*t+5*cos(10*t));
s3=hilbert(s3);
signal = generatesignal(1,1,10);
s0=signal(1:N);
S=[s0;s1;s2;s3];

[S0,Noise]=SNR(S(1:sk,:),sk,m,N,snr(1:sk));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=A*S0+Noise; 
% XI=A(:,2:sk)*S0(2:sk,:)+Noise; 
% RIN=XI*XI'/N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %����MVDR
a=array_line(0*pi/180,m);
R=X*X'/N;
invR=inv(R);
w1=invR*a/(a'*invR*a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ڶԽǼ��ص�MVDR��������10������������
RDL=R+10*eye(m);
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
pa=eye(m)-a*a'/(a'*a)+0.002*eye(m);
R3=pa*R*pa';
R3=(R3+R3')/2;
%% ����ȡ����ֵ��������
%����֪�źŸ���
[V,D] = eig(R3);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);
ed=D_sort(1:sk-1).*eye(sk-1);%������ֵ
us=V_sort(:,1:sk-1);%������ֵ��Ӧ����������
Aq=inv(pa)*us;
%% CVX���������
cvx_begin
variable w4(m,1) complex;
minimize((w4'*w4))
% minimize(norm(w4))
 subject to
  w4'*a==1;
  w4'*inv(pa)*us==zeros(1,sk-1);
cvx_end
%% ��ʽ��
% Aq=inv(pa)*us;
% B=inv(m*Aq'*Aq-Aq'*a*a'*Aq);
% w4=(a/m)*(a'*Aq*B*Aq'*a+1)-Aq*B*Aq'*a;

%%
sita=[-90:90];                    %% ɨ�跽��Χ
v=array_line(sita*pi/180,m);            %% ɨ�跽��ʸ��
B1=abs(w1'*v).^2; %MVER
B2=abs(w2'*v).^2; %
B3=abs(w3'*v).^2; %ESB
B4=abs(w4'*v).^2; %PROPOSED
%%
plot(sita,10*log10(B1/max(B1)),'k:','LineWidth',1.0);
hold on
plot(sita,10*log10(B2/max(B2)),'r--','LineWidth',1.0);
plot(sita,10*log10(B3/max(B3)),'g-.','LineWidth',1.0);
plot(sita,10*log10(B4/max(B4)),'b-','LineWidth',1.0);

hold on
% for i=1:sk  %����λ��
% plot([theta(i) theta(i)],[-90 0],'b:','LineWidth',1.0);
% end
xlabel('���﷽��(degree)');ylabel('����ͼ(dB)');
legend('MVDR','DL','ESB','Proposed algorithm');
grid on
% axis([-90 90 -90 0]);
% hold off


% pow=sum(abs(xs(1,:)).^2)/length(xs(1,:))
% 10*log10(pow)