function a=varray_line(doa,M,f,d)
%doa是各信号的方向组成的向量；M是阵元个数；f为载频，d为间距
c=3e8; ac=zeros(M,1);
    for k=1:M
        ac(k,1)=exp(-j*2*pi*f*d/c*sin(doa)*(k-1));
    end
    a=ac;
