function [S,Noise]=SNR(S,n,m,N,snr)
%S��Դ�ź�,n��Դ�źŸ�����m�����ߵ�Ԫ������N�ǲ���������snr�Ǹ��ź�����ȹ��ɵ�����(��λ��dB)��
Noise=randn(m,N)+j*randn(m,N);
for is=1:n
    S(is,:)=(S(is,:)-mean(S(is,:)))/std(S(is,:));%���ʹ�һ��
    S(is,:)=10^(snr(is)/20)*S(is,:);
end
for is=1:m
    Noise(is,:)=(Noise(is,:)-mean(Noise(is,:)))/std(Noise(is,:));
end