function [S,Noise]=SNR(S,n,m,N,snr)
%S是源信号,n是源信号个数，m是天线单元个数，N是采样点数，snr是各信号信噪比构成的向量(单位是dB)。
Noise=randn(m,N)+j*randn(m,N);
for is=1:n
    S(is,:)=(S(is,:)-mean(S(is,:)))/std(S(is,:));%功率归一化
    S(is,:)=10^(snr(is)/20)*S(is,:);
end
for is=1:m
    Noise(is,:)=(Noise(is,:)-mean(Noise(is,:)))/std(Noise(is,:));
end