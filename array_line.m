function A=array_line(a,m)
%a�Ǹ��źŵķ�����ɵ�������m����Ԫ��������Ԫ���벨��
count=length(a);
for is=1:count
    for js=1:m
        A(js,is)=exp(-j*pi*sin(a(is))*(js-1));
    end
end