function [S,Noise]=SNR(S,n,m,N,snr)
%S是源信号,n是源信号个数，m是天线单元个数，N是采样点数，snr是各信号信噪比构成的向量(单位是dB)。
Noise=randn(m,N)+j*randn(m,N);
Noise=powernorm(Noise);%功率归一化
S=powernorm(S);%功率归一化
for is=1:n    
    S(is,:)=sqrt(10^(snr(is)/10))*S(is,:);
end
