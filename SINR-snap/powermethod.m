function [m, u] = powermethod(A, ep, it_max)

%  求矩阵最大模特征值的幂法，调用格式为

%   [m, u] = powermethod(A, ep, it_max)

%  其中

%  A 为矩阵，ep 为精度要求，默认为1e-5,

%  it_max 为最大迭代次数，默认为100

%  m 为模最大的特征值，u 为 m 对应的特征向量

if nargin < 3 it_max = 100; end

if nargin <2 ep = 1e-5; end;

n = length(A); 
u = ones(n,1); %初始向量
k = 0; 
m1 = 0;

while k <= it_max

   v = A*u; 
   [vmax, i ] = max(abs(v));

   m  = v(i); %特征值
   u = v/m;%特征向量

   if abs(m-m1)<ep 
       break; 
   end

   m1 = m; 
   k = k+1;
end
