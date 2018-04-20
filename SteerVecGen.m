%--  102  --%
function [SteerVec] = SteerVecGen(TypeArray,Sita,Gama,M);
% 针对不同阵列生成导向矢量。
% 1、TypeArray == 'ULA'，生成均匀线阵(ULA)的导向矢量。
%                     (法线正向)
%                         ↑    ／
%                         ↑   ／
%                         ↑  ／    
%                         ↑ ／  θ
%                         ↑／   
%       M   M-1  ...  2   1  --->轴线正向
%    以阵元1为参考零点，d/λ固定为0.5。
%    以阵元M到阵元1为轴线正向，其逆时针转动角度为Sita，取值范围为[0，pi]。
%    参数要求：
%    (1) Sita是介于[0，pi]的方位角弧度；
%    (2) M是ULA的阵元数。
% 2、TypeArray == 'NUCA'，生成加芯均匀圆阵(NUCA)的导向矢量。
%    以下为一个5元NUCA阵列的排列示意图(俯视图)：
%                 3          /
%                           /θ
%            4    1    2     --->X轴正向 
%                       
%                 5   
%    阵元1位于圆心，并作为参考零点，R/λ固定为0.5(R为圆半径)。
%    参数要求：
%    (1) Sita是介于[0，  2pi]的方位角弧度；
%    (2) Gama是介于[0, 0.5pi]的俯仰角弧度。
%    (3) M是NUCA的阵元数。
% 3、TypeArray == 'UCA'，生成均匀圆阵(UCA)的导向矢量。
%    以下为一个4元UCA阵列的排列示意图(俯视图)：
%                 2          /
%                           /θ
%            3         1     --->X轴正向
%                       
%                 4   
%    M个阵元均匀排列在半径为R的圆圈上，以圆心-->阵元1为X轴正向，R/λ固定为0.5。
%    参数要求：
%    (1) Sita是介于[0，  2pi]的方位角弧度；
%    (2) Gama是介于[0, 0.5pi]的俯仰角弧度。
%    (3) M是UCA的阵元数。
% 4、TypeArray == 'UCA_1'，生成均匀圆阵(UCA)的导向矢量，以阵元1为参考零点。
%    以下为一个4元UCA阵列的排列示意图(俯视图)：
%                 2       /
%                        /θ
%            3         1     --->X轴正向
%                       
%                 4   
%    M个阵元均匀排列在半径为R的圆圈上，以圆心-->阵元1为X轴正向，R/λ固定为0.5。
%    参数要求：
%    (1) Sita是介于[0，  2pi]的方位角弧度；
%    (2) Gama是介于[0, 0.5pi]的俯仰角弧度。
%    (3) M是UCA的阵元数。
SteerVec = zeros([M,1]);
switch TypeArray
    case {'ULA'},
        mVec = (1:1:M).';   % 注意，添加了矢量转置操作，以保证输出是列矢量。
        SteerVec = exp(-j*(mVec-1)*pi*sin(Sita));  
    case {'NUCA'}
        SteerVec(1) = 1;
        mVec = (2:1:M).';   % 注意，添加了矢量转置操作，以保证输出是列矢量。
        SteerVec(2:M) = exp(-j*pi*sin(Gama)*sin((mVec-2)*2*pi/(M-1)-Sita));
    case {'UCA'}
        mVec = (1:1:M).';   % 注意，添加了矢量转置操作，以保证输出是列矢量。
        SteerVec = exp(-j*pi*sin(Gama)*sin((mVec-1)*2*pi/M-Sita));
    case {'UCA_1'}
        mVec = (1:1:M).';   % 注意，添加了矢量转置操作，以保证输出是列矢量。
        SteerVec = exp(-j*pi*sin(Gama)*sin((mVec-1)*2*pi/M-Sita));
        SteerVec = SteerVec/SteerVec(1);
    case {'UCA_2'}
        mVec = (1:1:M).';   % 注意，添加了矢量转置操作，以保证输出是列矢量。
        %///
        d_ratio = (9.52*sqrt(2)/2)/9.52;
        SteerVec = exp(-j*pi*d_ratio*sin(Gama)*sin((mVec-1)*2*pi/M-Sita));
        SteerVec = SteerVec/SteerVec(1);
end;
%-------------
