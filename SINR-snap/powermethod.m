function [m, u] = powermethod(A, ep, it_max)

%  ��������ģ����ֵ���ݷ������ø�ʽΪ

%   [m, u] = powermethod(A, ep, it_max)

%  ����

%  A Ϊ����ep Ϊ����Ҫ��Ĭ��Ϊ1e-5,

%  it_max Ϊ������������Ĭ��Ϊ100

%  m Ϊģ��������ֵ��u Ϊ m ��Ӧ����������

if nargin < 3 it_max = 100; end

if nargin <2 ep = 1e-5; end;

n = length(A); 
u = ones(n,1); %��ʼ����
k = 0; 
m1 = 0;

while k <= it_max

   v = A*u; 
   [vmax, i ] = max(abs(v));

   m  = v(i); %����ֵ
   u = v/m;%��������

   if abs(m-m1)<ep 
       break; 
   end

   m1 = m; 
   k = k+1;
end
