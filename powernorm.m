function [s1] = powernorm(s)
%s:������ź�
%s1�����ʹ�һ������ź�
[n,k]=size(s);
s1=zeros(n,k);
for i=1:n
pow = sum(abs(s(i,:)).^2)/length(s(i,:));
s1(i,:)=s(i,:)/sqrt(pow);%���źŹ��ʽ��й�һ��
end
end