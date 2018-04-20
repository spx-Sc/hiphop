function [S,Noise]=SNR(S,n,m,N,snr)
%S��Դ�ź�,n��Դ�źŸ�����m�����ߵ�Ԫ������N�ǲ���������snr�Ǹ��ź�����ȹ��ɵ�����(��λ��dB)��
Noise=randn(m,N)+j*randn(m,N);
Noise=powernorm(Noise);%���ʹ�һ��
S=powernorm(S);%���ʹ�һ��
for is=1:n    
    S(is,:)=sqrt(10^(snr(is)/10))*S(is,:);
end
