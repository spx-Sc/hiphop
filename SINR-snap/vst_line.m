function sta=vst_line(M,m,f,theta,fs,d)
%�˺�������������ʽ�������ʱ��ά��Ӧ������M����Ԫ������mʱ�����ӣ�
%d��Ԫ���
%thetaΪ�ź�����Ƕȣ���Ԫ���Ϊ�벨��,�����ź��ز�
a=varray_line(theta,M,f,d);
b=exp(j*2*pi/fs*f).^([0:m-1]');
sta=kron(a,b);