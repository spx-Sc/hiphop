function [s1] = powernorm(s)
%s:输入的信号
%s1：功率归一化后的信号
[n,k]=size(s);
s1=zeros(n,k);
for i=1:n
pow = sum(abs(s(i,:)).^2)/length(s(i,:));
s1(i,:)=s(i,:)/sqrt(pow);%对信号功率进行归一化
end
end