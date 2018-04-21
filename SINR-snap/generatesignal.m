function signal=generatesignal(Number_data,PRN,phaseshift)%数据个数是1，信号长度100ms
CAcode0=generateCAcode(PRN);


% phaseshift=round(1+1022*rand);
CAcode=[CAcode0(phaseshift:1023),CAcode0(1:phaseshift-1)];%相位延迟
% k=1:1023;
% figure(1)
% stairs(k,CAcode);
% axis([0 1024 -1.2 1.2]);
% grid on;

% ca_table1=[CAcode,CAcode,CAcode,CAcode];
% [a,b]=xcorr(ca_table1,'unbiased');
% figure(2)
% plot(b,a*1023);
% axis([-1100 1100 -70 1050])
% grid on;


%生成5位随机待发送二进制比特序列，并进行扩频编码
% Number_data=1;%待发送的二进制比特序列,修改数据位数N，则需要修改m(m序列的位数)
x_rand=rand(1,Number_data);
for i=1:Number_data
    if x_rand(i)>=0.5
        data(i)=1;  %data双极性数据码
    else data(i)=-1;
    end
end
% t_data=1:Number_data;
% figure(3)
% subplot(211)
% %stem(t,x)
% stairs(t_data,data);
% axis([0 7 -1.2 1.2])
% title('扩频前待发送二进制信息序列');
% subplot(212)
l=1:20460*Number_data;
y(l)=0;%y是20460*N个0
for i=1:Number_data
    for k=(20460*i-20459):(20460*i)
        y(k)=data(i);
    end
end
st_SpreadSpectrum_binary(l)=0;
for i=1:Number_data*20460
    if rem(i,1023)==0
        st_SpreadSpectrum_binary(i)=y(i)*CAcode(1023);
    else
        st_SpreadSpectrum_binary(i)=y(i)*CAcode(mod(i,1023));
    end
end
% t_SpreadSpectrum_binary=0:20460*Number_data-1;
% %stem(tt,s);
% stairs(t_SpreadSpectrum_binary,st_SpreadSpectrum_binary);
% axis([0 40920 -1.2 1.2]);
% title('扩频后的待发送序列码');

fc=9.46275e6;
fs=37.851e6;
fd=round(1000*(-10+20*rand));
dt=1/fs;
f_data=50;
T=Number_data/f_data;%100ms
t=0:dt:T-dt;
st_afterSS_binary=rectpulse(st_SpreadSpectrum_binary,length(t)/length(st_SpreadSpectrum_binary));%将冲激信号补成矩形信号
s_bpsk_afterSS=st_afterSS_binary.*cos(2*pi*(fc+fd)*t);
% figure(4)
% plot(ts,s_bpsk_afterSS);
% xlabel('s');
% axis([0.5/1023000 9.5/1023000 -1.2 1.2])
% %axis([0 10.5 -1.2 1.2])
% title('扩频后bpsk信号时域波形');

% figure(5)
% ff=(0:length(t)-1)*fs/length(t);
% y=abs(fft(s_bpsk_afterSS))*2/length(t);
% plot(ff(0,length(t)/2-1),y(0,length(t)/2-1));
% grid on;
% axis([7e6 12e6 -0.2 1])
%noise=round(-20+10*rand);%-20到-10dB的噪声
signal=s_bpsk_afterSS;
%signal=awgn(s_bpsk_afterSS,noise);
end



