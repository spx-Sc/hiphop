function a=varray_line(doa,M,f,d)
%doa�Ǹ��źŵķ�����ɵ�������M����Ԫ������fΪ��Ƶ��dΪ���
c=3e8; ac=zeros(M,1);
    for k=1:M
        ac(k,1)=exp(-j*2*pi*f*d/c*sin(doa)*(k-1));
    end
    a=ac;
