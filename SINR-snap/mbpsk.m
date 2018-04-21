%OOK,2PSK,文件名binarymod.m
function [st,ss]=mbpsk(fc,N_sample,N,B)
%fc载频，N_sample每个码元抽样个数，fs=Ts/fc/N_sample
%Ts=2/B，B为带宽，Ts为码元速率，N 为码元个数。
A=1;
% fc=100;   
% N_sample=6;
% N=50;%码元数
% Ts=1;
Ts=2/B;
dt=Ts/fc/N_sample;  %波形采样间隔
t=0:dt:N*Ts-dt;
Lt=length(t);
%产生二进制信源
d=sign(randn(1,N));
dd=sigexpand((d+1)/2,fc*N_sample);
gt=ones(1,fc*N_sample); %NRZ波形
d=sign(randn(1,N));
dd=sigexpand((d+1)/2,fc*N_sample);
gt=ones(1,fc*N_sample); %NRZ波形
d_NRZ=conv(dd,gt);
ht=A*cos(2*pi*fc*t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1);
%2PSK信号
d_2psk=2*d_NRZ-1;
s_2psk=d_2psk(1:Lt).*ht;
% plot(t,s_2psk);
% axis([0 10  -1.2 1.2]);
% xlabel('2PSK');
% s_2psk=hilbert(s_2psk);
% figure(2)
% N=length(t);
% fs=1/dt;
% ff=[0:N-1]*fs/N;
% plot(ff,20*log10(abs(fft(s_2psk))));
% grid on;
% axis([fc-4/Ts fc+4/Ts  0 60]);
% xlabel('2PSK频谱(dB/Hz)');
st=t;
ss=s_2psk;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

