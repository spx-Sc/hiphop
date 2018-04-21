function sta=vst_line(M,m,f,theta,fs,d)
%此函数产生线阵形式下虚拟空时二维响应，其中M是阵元个数，m时间因子，
%d阵元间距
%theta为信号入射角度，阵元间距为半波长,入射信号载波
a=varray_line(theta,M,f,d);
b=exp(j*2*pi/fs*f).^([0:m-1]');
sta=kron(a,b);