function A=array_line(a,m)
%a是各信号的方向组成的向量；m是阵元个数；阵元间距半波长
count=length(a);
for is=1:count
    for js=1:m
        A(js,is)=exp(-j*pi*sin(a(is))*(js-1));
    end
end