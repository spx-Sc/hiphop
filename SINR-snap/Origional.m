% function SINROutVec = ClassicalWeinerForKongpin(SigSitaVec,SigGamaVec,P)
%master bench   %696
clear all;
clc;
format long e
% 程序开始，显示程序开始时间。
disp(sprintf('------程序开始仿真 ......'));
%%%%%%%%%%%%%%%%%%%%%%%%%%

SigSitaVec = (pi/180)*[10];       % 信号来向
Jam=[50;10;-30];
JamSitaVec = Jam*pi/180; %干扰方位角  平面阵列
AlgorithmType = '5';                                      %选择算法执行类型
ybsn=0.0001;
C=1;
dirta=2;
P=120;

Num_Sig =size(SigSitaVec);                                                                 %GPS卫星信号数目
Num_Jam = size(Jam);                                                              %干扰个数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_begin = cputime;
Time_clock = fix(clock);
disp(sprintf('程序开始时间： %d-%d-%d  %d:%d:%d',...
              Time_clock(1),Time_clock(2),Time_clock(3), ...
              Time_clock(4),Time_clock(5),Time_clock(6)));
FlagDisp = 1;
%FlagFixLen = 1;
TypeArray = 'ULA';  %%%天线阵列形状选择，针对不同阵列生成导向矢量，详见SteerVecGen（）函数
%  SttrandBegin =round(rand*100000);
SttrandBegin = 70;
rand('state',SttrandBegin);
SttrandConfg = round(rand*10000);%产生一个4位数的随机数
SttrandBit   = round(rand*10000);
SttrandnJam  = round(rand*10000);   %干扰？    
SttrandnNoi  = round(rand*10000);   %噪声？

%JamSitaVec = [70;90;150;130]*pi/180; %干扰方位角   均匀线阵
% JamSitaVec = [-30;30;70]*pi/180; %干扰方位角  平面阵列
JamGamaVec = [45;45; 45]*pi/180;%干扰俯仰角
% SigSitaVec = 10;%信号方位
% SigGamaVec = 30;%信号俯仰

disp(sprintf(' 1、基本参数设定 '));
%% % 1、基本参数设定 --------------------------------------------%%%%%%
%% %    1.1 天线阵列和联合空时处理的时间抽头个数 --------------------%%%%
rand('state',SttrandConfg);            % 装定生成配置用伪随机数发生器的状态,如：执行rand('state',0);rand(1)   rand('state',SttrandConfg);rand(1)
M = 10;                  % 阵元个数
Ppie =1;                  %时域抽头个数
if ((Ppie-1)/2-round((Ppie-1)/2))>10^(-12),
    disp('▲告警▲要求时域抽头个数 P 是奇整数。');
else
    P_h = round((Ppie+1)/2);
    DlySTAP = P_h - 1;  % STAP带来的延迟,虽然可以处于[0，P-1]之间，但目前固定为2。%这不太懂为什么要这么做
end;

%%    1.2 采样频率、采样点数、抽取倍数 ---------------------------%%%
Rc = 1.023*10^6;        % GPS C/A码的码片速率
F0 = 1540*Rc;           % GPS信号的L1频率，射频频率
Fc =   70*Rc;           % 载波频率，采用带通采样方法
G  =   40*2;              % 每一个码片周期内的采样位数，要求采样频率必须是码片速率的偶整数倍
G_h = G/2;
G_q = G/4;
Q  =    8;              % 抽取倍数
R  = G / Q;             % 抽取后每个码片时间周期内的采样点数（快拍数）
R_h = R/2;
Fsam =  G*Rc;           % 要求采样频率是码片速率的偶整数倍
Num_Bit = 2;            % 信息比特的数目
Len_PRN = 1023;         % 每个扩频码长1023个码片
Len_PRNSnap = Len_PRN*R;% 每个完整的PRN码序列的采样快拍数
Len_PRNSam  = Len_PRN*G;% 每个完整的PRN码序列的采样点数
Num_Chip = Num_Bit  * Len_PRN;      % 码片长度
Num_Sam  = Num_Chip * G;            % 采样序列的长度
Num_Snap = Num_Chip * R;            % 联合空时处理后的快拍数，最终阵列处理输出
                                    % 的数据长度
%%    1.3 接收机天线处白噪声方差(功率) ---------------------------%%%
NoiPow = 1;                         % 热噪声的功率
%%    1.4 GPS信号的功率、幅度和初相 -----------------------------%%%

switch TypeArray                    % 定义各个卫星信号的方位角和俯仰角
    case {'ULA'}
%         SigSitaVec = (pi/180)*[0];% 各个卫星的方位角
        SigGamaVec = (pi/180)*[45];% 各个卫星的俯仰角，对于线阵，
                                                        % 该取值无意义
    case {'NUCA','UCA','UCA_1','UCA_2'}
        SigSitaVec = (pi/180)*[ 150; 170;30;200;251;307];% 各个卫星的方位角
        SigGamaVec = (pi/180)*[ 45; 45; 45; 45; 45; 45];% 各个卫星的俯仰角
end;
SVIDVec = [1;5;9;13;17;19];               % 各个卫星序号
SVTaoRatioVec = [20+round(50*rand([Num_Sig,1]))];   % 各个卫星相对于理想零点的延迟/Tsnap。    round()四舍五入
SNRVec = -27.5 * ones([Num_Sig,1]); % 各个卫星信号的信噪比(dB值)，取为GPS   ones([Num_Sig,1])信号个数行，有几个信号就有几个信噪比
                                    % C/A码的典型信噪比，这里多减7dB是为了适合
                                    % 后面的STAP前SNR数值达到-21.9dB。                                    
SigPowVec = NoiPow*10;% 各个GPS信号的接收功率。由信噪比和噪声功率，反推信号功率

SigFaiVec = rand([Num_Sig,1])*2*pi; % 有用信号载波的初相，由于接收机正交数字载
                                    % 波的初相是固定的，因此可以认为发射信号的    不太明白？？？？？
                                    % 载波初相是随机分布的。
                                   
%%    1.5 干扰的数目、类型、功率和入射角 -------------------------%%%
 JamKindVec = ['FBJam';'FBJam';'FBJam'];    % 干扰类型
%JamKindVec = ['FBJam';'FBJam';'FBJam';'STJam'];    % 干扰类型
                                       % ='FBJam' 宽带噪声干扰；
                                       % ='NBJam' 窄带噪声干扰；
                                       % ='STJam' 单频干扰；
JNRVec = [30;30;30];                % 干噪比，干扰相对于噪声的功率比(dB值)
%JNRVec = [40;40;35;35;35]-7.5;        % 干噪比，干扰相对于噪声的功率比(dB值)
%JamPowVec = NoiPow*10.^(JNRVec/10);    % 各个干扰分量的功率(这里指的是干扰的复功
JamPowVec = NoiPow*10.^(JNRVec/10);          % 率,实功率是复功率的一半)  
                                       
JamFrqOffVec = [0; -1; 0.5;-1;+1]*Rc;          
                                       % 干扰中心频率相对于载波的频率差,仅对音频
                                       % 干扰有意义。其定义为(fi-Fc)，Fc为信号载
                                       % 波，对于基带生成情况，Fc=0。
                                       % 考虑到中频带通滤波器的通带范围，因此干扰
                                       % 频差系数的绝对值不能超过1。
JamFaiVec = rand([Num_Jam,1])*2*pi;    % 干扰的载波初相

%%    1.6 带通滤波器,它是一个通带中心频率在0.25*Fsam的带通滤波器 ---%%%
Rp = 1;                             % 通带最大允许纹波
Rs = 80;                            % 阻带衰减
% f1 = [0.2/15.6  6.5/15.6  9.1/15.6  15.4/15.6];    % 通、阻带截止频率设定
f1 = [0.1/G_h  (G_q-1)/G_h  (G_q+1)/G_h  (G_h-0.1)/G_h];
                                    % 通、阻带截止频率设定
a1 = [    0    1    0     ];        % 指定通阻带特性
dev1 = [10^(-Rs/20)  (10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)]; 
                                    % 根据通带纹波和阻带衰减准备参数
[n,fo,ao,wo] = firpmord(f1,a1,dev1); % 确定FIR滤波器的阶数
n1 = ceil(n/2)*2 + 10;      % 采用firpm.m设计FIR滤波器往往会稍微达不到阻带
                            % 衰减要求，故提前多加一些长度。而且，要求n1是
                            % 偶数，以便后面的群延迟n1/2是一个整数。
BPFVec = firpm(n1,fo,ao,wo); % 产生线性相位的带通滤波器
GrpDlyBPF = round((n1-1)/2);    % 线性带通滤波器的群延迟
%%%%--------------------------------------------------------------------%%%
%%    1.7 低通滤波器设计 ---------------------------------------%%%
%            用于正交下变频时I、Q两路混频后数据的低通滤波，并准备抽取
Rp = 1;                             % 通带最大允许纹波
Rs = 50;                            % 阻带衰减
% f2 = [1.1/15.6  2.5/15.6];           % 通、阻带截止频率设定
f2 = [1/G_h  R_h/G_h];                % 通、阻带截止频率设定
a2 = [ 1    0  ];                   % 指定通阻带特性
dev2 = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)]; 
                                    % 根据通带纹波和阻带衰减准备参数
[n,fo,ao,wo] = firpmord(f2,a2,dev2); % 确定FIR滤波器的阶数
n2 = ceil(n/2)*2 + 4;%10;      % 采用firpm.m设计FIR滤波器往往会稍微达不到阻带
                            % 衰减要求，故提前多加一些长度。而且，要求n2是
                            % 偶数，以便后面的群延迟n2/2是一个整数。
LPFVec = firpm(n2,fo,ao,wo);         % 产生线性相位的低通滤波器
GrpDlyLPF = round((n2-1)/2);        % 线性低通滤波器的群延迟

FlagLoJit = 0;                      % 仿真时考虑时序抖动问题
% LoJit_std = 2*10^(-12);             % 本振频率和ADC采样频率的合成等效
%                                     % 时序抖动的方差
% if FlagLoJit == 1,
%     disp(sprintf('  (9)仿真中考虑了本振频率和ADC采样频率的时序抖动，二者合成等效的时序抖动的方差为 %4.2g ps。',...
%                    LoJit_std*10^12) );
% end;
%%%%-------------------------------------------------------------------%%%%
%%   2、仿真数据准备 ------------------------------------------%%%%
 disp(sprintf(' 2、从中频开始生成数据，经混频和低通滤波、抽取后得到基带数据。'));
%%     2.1 基本配置情况说明 ------------------------------------%%%
        Num_BitRaw = Num_Bit + ceil((GrpDlyBPF+GrpDlyLPF)/Len_PRNSam/2)*2;
                                                % 原始信息比特个数，要求是偶数
        Num_ChipRaw = Num_BitRaw * Len_PRN;     % 原始码片个数
        Num_SamRaw  = Num_ChipRaw * G;          % 原始采样个数
        Num_SamADC = Num_SamRaw - GrpDlyBPF;    % ADC采样得到的数据序列长度
        nVec = ( 1:1:Num_SamRaw ).';            % 采样序号列矢量，从1开始 
%%     2.2 针对C/A码情况生成码片数据和码片采样数据 ----------------%%%
        Ratio_Sig  = Fc/Fsam;                   %一个码片采样点，需要几个载波周期
        BitMat_Raw = zeros([Num_BitRaw,Num_Sig]);
        PRNMat     = zeros([ Len_PRN  ,Num_Sig]);
        DSSnapMat  = zeros([round(Num_SamRaw/Q) ,Num_Sig]);
        SigVecMat  = zeros([Num_SamADC,Num_Sig]);
        SigIFallMat= zeros([Num_SamADC,   M   ,Num_Sig]);%M = 4; % 阵元个数 Num_Sig 信号个数
        SigIFMat   = zeros([Num_SamADC,   M   ]);
        SigSteerallMat= zeros([M,Num_Sig]);
        for k=1:Num_Sig,
            rand('state',SttrandBit+(k-1)*237); % 均匀分布伪随机数状态装定
            BitMat_Raw(:,k) = [ AJ_AryBalance(rand([           Num_Bit,1]));  ...
                                AJ_AryBalance(rand([Num_BitRaw-Num_Bit,1])) ] ;
                                    % 产生一个+1，-1个数相等的均衡序列，而且
                                    % 要求前NBit个序列元素之和也等于0。
            [G1Vec,PRNMat(:,k)] = AJ_GPSCAPRNGen(SVIDVec(k)); 
                                    % 给出第k颗GPS卫星的PRN码，PRNMat为长度1023码片的CA码
            DSChipVec = kron( BitMat_Raw(:,k),PRNMat(:,k) );
            DSSnapVec_NoTao = kron(DSChipVec,ones([R,1]));
                                    % 对GPS信号码片的快拍采样序列矢量，尚未
                                    % 叠加延迟
            DSSnapMat(:,k) = circshift(DSSnapVec_NoTao,SVTaoRatioVec(k));
                                    % 把各个卫星PRN码的相对延迟体现出来
            DSSamVec = kron(DSSnapMat(:,k),ones([Q,1]));
            tempVec1 = exp( j*(2*pi*Ratio_Sig*(nVec-1)+SigFaiVec(k) ...
                                     +0.5*(1- DSSamVec)*pi)        ) ;
                                    % 码片+1对应于叠加相位0，码片-1对应于叠加
                                    % 相位pi。功率稍后指定。
                                    %tempVec1为模拟下变频后的中频信号，%其中SigFaiVec信号初相，此时已经完成了中频调制
            SigVecMat(:,k) = IFFltandPowSet(tempVec1,SigPowVec(k),BPFVec);
                                    %对模拟下变频后的中频信号进行模拟带通滤波，其中SigPowVec信号功率
            SigSteerVec = SteerVecGen(TypeArray,SigSitaVec(k),SigGamaVec(k),M);
            Alb=SigSteerVec;
                                    %针对不同阵列生成导向矢量，TypeArray阵元类型，信号的方位俯仰，阵元个数
            SigSteerallMat(:,k) = SigSteerVec;
                                    %生成信号个数Num_Sig个阵列流型
            SigIFallMat(:,:,k) = kron(SigVecMat(:,k),SigSteerVec.');
                                    %阵列导向矢量与中频信号相乘
            SigIFMat = SigIFMat + SigIFallMat(:,:,k);
                                    % 得到每个天线对有用信号的复数接收数据,并
                                    % 将其累加起来，作为有用信号
        end;
        clear tempVec1
        clear SigVecMat
        
%%     2.3 干扰生成 -------------------------------------------%%%
        randn('state',SttrandnJam);             % 高斯分布伪随机数状态装定
        JamIFMat = zeros([Num_SamADC,M]);       % 天线阵接收到的干扰复数数据
        JamSteerallMat = zeros([M,Num_Jam]);
        for k=1:Num_Jam,
            switch JamKindVec(k,:),
                case {'FBJam'},                 % 干扰是高斯白噪声序列
                    tempVec1 =    randn(Num_SamRaw+GrpDlyLPF,1) ...
                               +j*randn(Num_SamRaw+GrpDlyLPF,1) ;
                                                % 生成一个复数白噪声序列
                    tempVec0 = filter(LPFVec,1,tempVec1);
                    Ratio_Jam = Ratio_Sig;      % 满带噪声干扰的中心频率与信号相同。
                    tempVec2 =    tempVec0((GrpDlyLPF+1):(GrpDlyLPF+Num_SamRaw)) ... 
                              .* exp(j*(2*pi*Ratio_Jam*(nVec-1)+JamFaiVec(k)));
                case {'PBJam'},
                  
                case {'STJam'},                 % 干扰是单频干扰
                    Ratio_Jam = (Fc+JamFrqOffVec(k))/Fsam;
                    tempVec2 = exp(j*(2*pi*Ratio_Jam*(nVec-1)+JamFaiVec(k)));
                                                % 功率随后会赋予
            end;
            JamVec = IFFltandPowSet(tempVec2,JamPowVec(k),BPFVec);
            SteerVec_Jam = SteerVecGen(TypeArray,JamSitaVec(k),JamGamaVec(k),M);
                                            % 各个天线接收到的信号附加相位差
                                            % 表现为对接收信号的模1复数相乘。
%             JamSteerallMat(:,k) = SteerVec_Jam;
            JamIFMat = JamIFMat + kron(JamVec,SteerVec_Jam.');
                                            % 将所有宽带干扰的复数接收数据叠加起来        
        end;
        clear tempVec2
%%     2.4 白噪声生成 -----------------------------------------%%%
        randn('state',SttrandnNoi);         % 高斯分布伪随机数状态装定
        NoiIFMat = sqrt(NoiPow/2)*(      randn(Num_SamADC,M) ...
                                     + j*randn(Num_SamADC,M) );
%%    2.5 信号生成情况说明 ------------------------------------%%%
        
%%    2.6 信号、干扰和白噪声合成为ADC输出的数据  -----------------%%%
        RecIFMat = SigIFMat + JamIFMat + NoiIFMat;        
                                            % 得到天线阵的复数接收数据矩阵   
%         if FlagLoJit == 1,                  % 判断仿真时是否考虑时序抖动问题
%             LoJitMat = exp(j*2*pi*F0*LoJit_std*randn(size(RecIFMat)));
%             RecIFMat = RecIFMat.*LoJitMat;
%         end;
      
        RecADCMat = real(RecIFMat);         % 实际中ADC采样得到的是实数序列        
        SigADCMat = real(SigIFMat);         %目标信号实数数据
        JamADCMat = real(JamIFMat);         %干扰实数数据
        NoiADCMat = real(NoiIFMat);         %噪声实数数据
        SigADCallMat = real(SigIFallMat);   %信号个数没有叠加的数据
        
%%    2.7 混频和低通滤波(需要按照抽取的要求确定阻带截止频率)-------%%%
        Fai_Mix = 0;                    % 本地混频载波初相
        RecBFMat = MixandLPF(RecADCMat,Ratio_Sig,Fai_Mix,LPFVec);%% 对ADC采样得到的带通数据序列进行混频和低通滤波，并且将I、Q两路输出信号重新组合为一个复数的基带信号序列。
        SigBFMat = MixandLPF(SigADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        JamBFMat = MixandLPF(JamADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        NoiBFMat = MixandLPF(NoiADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        SigBFallMat = zeros([Num_SamADC-GrpDlyLPF,M,Num_Sig]);
        for k=1:Num_Sig,
            SigBFallMat(:,:,k) = MixandLPF(SigADCallMat(:,:,k),Ratio_Sig,Fai_Mix,LPFVec);
        end;
%%     2.8 抽取 -----------------------------------------------%%%
        nVec3 = ( (0:1:Num_Snap-1)*Q+1 ).'; % Q倍抽取的位置号构成一个矢量
        RecSnapMat = RecBFMat(nVec3,:);     % 接收信号混频滤波后的Q倍抽取
        SigSnapMat = SigBFMat(nVec3,:);     % 有用信号混频滤波后的Q倍抽取
        JamSnapMat = JamBFMat(nVec3,:);     % 干扰混频滤波后的Q倍抽取
        NoiSnapMat = NoiBFMat(nVec3,:);     % 噪声混频滤波后的Q倍抽取
        SigSnapallMat = SigBFallMat(nVec3,:,:);
%% %----------------------仿真信号生成完毕-----------------------%%%
%%  3、自适应滤波 ---------------------------------------------%%%%
%%    3.1 基本参数设定和数据准备 --------------------------------%%%
disp(' 3、开始求解权值矢量');
RecSTAPMat = zeros([M*Ppie,Num_Snap-(Ppie-1)]);
% 以x(m,n)表示第m个天线支路上第n个时刻的快拍值，则第n个时刻STAP滤波列矢量规定为：
% [x(1,n+P-1); ... x(M,n+P-1);x(1,n+P-2); ... x(M,n+P-2); ...;x(1,n); ...
% x(M,n)]       
for n=1:(Num_Snap-(Ppie-1)),
    RecSTAPMat(:,n) = reshape(RecSnapMat((n+Ppie-1):-1:n,:).',M*Ppie,1); %真实的接收信号数据
    SigSTAPMat(:,n) = reshape(SigSnapMat((n+Ppie-1):-1:n,:).',M*Ppie,1); %第一种需要用的量，目标信号数据   
end;

%%    3.2 滤波权值计算 ----------------------------------------%%%
% AlgorithmType = '5';                                      %选择算法执行类型
%'0'：功率倒置算法
%'1'：MVDR算法
%'2'：经典维纳滤波算法
%'3'：CSA_MSNWF:相关相减结构的多级维纳滤波    P72
%'4'：多步迭代降秩算法    P73
disp(sprintf('----执行第%s种算法 ',AlgorithmType));
Time_begin1 = cputime;
%%%%%%%%%%%%%专门用于计算Rx%%%%%%%%%%%%%%%%%%
        L =1024;      
        xAnaMat = RecSTAPMat(:,1:L);%输入信号
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %输入向量得自相关矩阵
        Rx=RxMat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch AlgorithmType, 
      %% 功率倒置算法  
        case {'0'},
        disp('    ——采用功率倒置算法进行权值计算,');
        L =1024;
        xAnaMat = RecSTAPMat(:,1:L);%输入信号
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %输入向量得自相关矩阵
        Rx=RxMat;
        a=zeros([M*Ppie,1]);
        a(1)=1;  % 功率倒置
        wVec = (inv(RxMat)*a)/(a'*inv(RxMat)*a);
  %%  MVDR算法
        case {'1'},
        disp('    ——采用MVDR算法进行权值计算,');
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);%输入信号
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %输入向量得自相关矩阵
        a=zeros([M*Ppie,Num_Sig]);
        %a=st_cir(SigSteerVec,P,Fsam,Fc); %用一颗星导向矢量
        for k=1:Num_Sig
             a(:,k)=st_cir(SigSteerallMat(:,k),Ppie,Fsam,Fc);
        end     
        I_sig=ones(Num_Sig,1);
        wVec = (inv(RxMat)*a)/(a'*inv(RxMat)*a)*I_sig;
        %% %-------------------第1种经典维纳滤波-------------------------%%
    case {'2'},
        disp('    ——采用经典维纳滤波算法进行权值迭代,');
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);%输入信号
        [MM,NumWeiner] = size(xAnaMat);
        dRow    = SigSTAPMat( 1 ,1:L);%期望信号
      
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %输入向量得自相关矩阵
        pxVec = xAnaMat*dRow'   /NumWeiner;  %期望信号与输入信号的互相关矢量
        wVec = zeros([MM,1]);
        wVec =   (RxMat \ pxVec); % 采用矩阵左除，通过高斯消元法求解维纳最优解。左除式A\B，则相当于inv（A）*B
        wVec = wVec/norm(wVec);
%% %----------------------%第2种CSA_MSNWF_wVec:相关相减结构的多级维纳滤波器%----%
    case {'3'},
%         TypeConst = 'SimpConst';
        disp('    ——采用CSA_MSWNF算法进行权值迭代,');
        L = 1024;    % 用L个个XX'求平均，估计hi和wi，最大
                                        % 不超过128个数据点。要求L是2的整数
                                        % 次幂，否则后面求平均值的操作会带来
                                        % 误差。
        xAnaMat = RecSTAPMat(:,1:L);
        TypeConst = 'SSSCConst';
        wVec = CSA_MSNWF_wVec(M,Ppie,L,TypeConst,SigSteerVec,xAnaMat);                               
                                        

%% %-----------------------------第3种：ICSA_MSNWF_wVec:多步迭代相关相减降秩算法--------------------------%%
    case {'4'},
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);
        TypeConst = 'SSSCConst';
        wVec = ICSA_MSNWF_wVec(M,Ppie,L,TypeConst,SigSteerallMat(:,1),xAnaMat);
    case{'5'}
        Kmax=9999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
        g=1; % 我们期望的反应
%         sitai(1:40)=linspace(-90,90,40);
%         sitai(41:80)=linspace(-90,-5,40);
%         sitai(81:120)=linspace(5,90,40);
%         for ff=1:120
%         sitai(1:20)=linspace(-90,90,20);
%         sitai(21:40)=linspace(-90,-5,20);
%         sitai(41:60)=linspace(5,90,20);
        sitai=linspace(-90,90,60);
        for ff=1:60
            SteerVec_sitai(ff,:) = SteerVecGen(TypeArray,sitai(ff),45*pi/180,M);
        end
        as= SteerVec_sitai;
        
%         for i=1:60
%             aba(2*i-1,1)=real(as(i,:));  %aba is the a-ba; aba(120,2)
%             aba(2*i,1)=imag(as(i,:));
%             aba(2*i+119,1)=imag(as(i,:));
%             aba(2*i+120,1)=-real(as(i,:));
%         end
        for i=1:60
            denglou(:,i)=[real(as(i,:)') ; imag(as(i,:)')];
            denglou(:,i+60)=[imag(as(i,:)') ; -real(as(i,:)')];
        end
% 　　i=1;
%         aba1=[real(as(i,:)') ; imag(as(i,:)')];
%         aba2=[imag(as(i,:)') ; -real(as(i,:)')];
%         for i=2:60^^
%             aba1=[aba1 ; real(as(i,:)') ; imag(as(i,:)')];
%             aba2=[aba2 ; imag(as(i,:)') ; -real(as(i,:)')];
%         end
%         denglou=[aba1;aba2];
      

Rxp=[real(Rx),-imag(Rx);imag(Rx),real(Rx)];   %Rxp!=Rx!!!!!!!!!!!!!
         
        
        
%         for i=1:120
%             if abs(sitai(i))>dirta
%                 di(i)=0;
%             else
%                 di=g;
%             end
%         end
        for i=1:60
            if abs(sitai(i))>dirta
                di(i)=0;
            else
                di(i)=g;
            end
            di(i)=real(di(i));
            if abs(sitai(i))>dirta
                di(i+60)=0;
            else
                di(i+60)=g;
            end
            di(i+60)=imag(di(i+60));
        end
       % Rx=RxMat;
        
        
%         u(1,:)=(di);   %u's tructure : u(k,i) ,it should be in the outer of loop
%         for K=1:Kmax
%             w(:,K)=0;
%             for i=1:120
%                 if u(K,i)<ybsn
%                     fi(i)=0;
%                 else
%                     fi(i)=C/u(K,i);
%                 end
%                 Df(i,i)=fi(i);
%             end
u=(di);   %u's tructure : u(k,i) ,it should be in the outer of loop
for K=1:Kmax
    w(:,K)=zeros(20,1);
    for i=1:120                       %using method m=1 another place used
        if u(k,i)<ybsn
            fi(i)=0;
        else
           % fi(i)=C/u(k,i);   j=1
           fi(i)=2*C*(u(k,i)-ybsn)/u(k,i);    %m=2version
        end
        Df(i,i)=fi(i);
    end
    %%%%%
    ws=(Rxp+denglou*Df*denglou')\(denglou*Df*di');%%%the problem still exist.continue!
    yita(K)=1;
    %%%
    
    w(:,K+1)=w(:,K)+yita(K)*(ws-w(:,K));
    for i=1:120    % warning!
        aba=denglou(:,i);
        u(K+1,i)=abs(di(i)-w(:,K+1)'*aba);
    end
    while (  jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)>=jhs(w(:,K),Rxp,u(K,:),bi,ybsn,C)   )
        yita(K)=yita(K)-0.1;                 % how fast the number is decreadsed in not sure
        w(:,K+1)=w(:,K)+yita(K)*(ws-w(:,K));
        for i=1:120    % warning!
            aba=denglou(:,i);
            u(K+1,i)=abs(di(i)-w(:,K+1)'*aba);
        end
          jhs(w(:,K+1),Rxp,u(K+1,:),bi,ybsn,C)
    end
    % form now, the value of k is virtually increased.
    kjia=K+1;
    for i=1:120    % warning!
        aba=denglou(:,i);
        u(kjia,i)=abs(di(i)-w(:,kjia)'*aba);
    end
    %%%%%%%%
    for i=1:120     %using method m=1   another place used
        if u(kjia,i)<ybsn
            fi(i)=0;
        else
            fi(i)=C/u(kjia,i);
        end
    end
    %%%%%%%%%
%     jhs(w(kjia+1),Rxp,u(kjia+1,:),bi,ybsn,C)
%     jhs(w(:,kjia),Rxp,u(kjia,:),bi,ybsn,C)
   % if ((jhs(w(kjia+1))-jhs(w(:,kjia)))<minchaju )
    if ( jhs(w(:,kjia),Rxp,u(kjia,:),bi,ybsn,C)- jhs(w(:,kjia-1),Rxp,u(kjia-1,:),bi,ybsn,C))<minchaju;
        break;
    end
end
for i=1:10
    wVec(i)=w(i,K+1)+w(i+10,K+1)*j;
end
wVec=wVec';
    %% %-----------------------------------------------------------%%
end      %算法类型选择结束点


        wVecallMat = kron(ones([1,Num_Sig]),wVec);  
%% %--------------------------------------------------------%%
Time_mid = cputime;
Time_used = Time_mid - Time_begin1;    
disp(sprintf('----求权所用时间： Time_used = %g 秒。',Time_used));
%% %%%% 4、滤波,求出阵列输出信号 ------------------------------------------%%%%%%
%% %    4.1 滤波 -----------------------------------------------%%%%
disp(' 4、权值矢量计算完毕，现在开始用该权值矢量进行抗干扰滤波');
yOutallMat   = zeros([Num_Snap,Num_Sig]);%NumSnap抽取后的数据点数
JamOutallMat = zeros([Num_Snap,Num_Sig]);%干扰
NoiOutallMat = zeros([Num_Snap,Num_Sig]);%噪声
SigOutallMat = zeros([Num_Snap,Num_Sig]);%有用信号
for k=1:Num_Sig
    yOutallMat(:,k)   = STAPFilter(RecSnapMat,Ppie,wVec,DlySTAP);%真实数据，DlySTAP：stap延迟
    SigOutallMat(:,k) = STAPFilter(SigSnapallMat(:,:,k),Ppie,wVec,DlySTAP);%有用信号
    JamOutallMat(:,k) = STAPFilter(JamSnapMat,          Ppie,wVec,DlySTAP);%干扰
    NoiOutallMat(:,k) = STAPFilter(NoiSnapMat,          Ppie,wVec,DlySTAP);%噪声
end;
%% %--------------------------------------------------------%%
Time_mid = cputime;
Time_used = Time_mid - Time_begin1;    
disp(sprintf('----滤波所用时间： Time_used = %g 秒。',Time_used));
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %    4.2 计算STAP输入和STAP输出SINR ---------------------------%%%%
%%% 计算理想情况下，针对每一个天线的最佳权值，所得到的最优输出信干噪比
%%%%%%%%%%%%%%%%%输入%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowSigInVec    = zeros([Num_Sig,1]);
PowJamInVec    = zeros([Num_Sig,1]);
PowNoiInVec    = zeros([Num_Sig,1]);
PowJamNoiInVec = zeros([Num_Sig,1]);
for k=1:Num_Sig,
    PowSigInVec(k)    = mean(abs(               SigSnapallMat(:,1,k)).^2);%输入信号功率
    PowJamInVec(k)    = mean(abs(               JamSnapMat(:,1)    ).^2);%输入干扰功率
    PowNoiInVec(k)    = mean(abs(               NoiSnapMat(:,1)    ).^2);%输入噪声功率
    PowJamNoiInVec(k) = mean(abs(RecSnapMat(:,1)-SigSnapallMat(:,1,k)).^2);%输入干扰加噪声功率

end;
SNRInVec     = 10*log10(PowSigInVec./PowNoiInVec);%输入信噪比
SINRInVec    = 10*log10(PowSigInVec./PowJamNoiInVec);%%输入信干噪比
JSRInVec     = 10*log10(PowJamInVec./PowSigInVec);%输入干信比
%%%%%%%%%%%%%%%%输出%%%%%%%%%%%%%%%%%%%%%%%%
PowSigOutVec    = zeros([Num_Sig,1]);
PowJamOutVec    = zeros([Num_Sig,1]);
PowNoiOutVec    = zeros([Num_Sig,1]);
PowJamNoiOutVec = zeros([Num_Sig,1]);
for k=1:Num_Sig,
    PowSigOutVec(k)    = mean(abs(                SigOutallMat(:,k)).^2);
    PowJamOutVec(k)    = mean(abs(                JamOutallMat(:,k)).^2);
    PowNoiOutVec(k)    = mean(abs(                NoiOutallMat(:,k)).^2);
    PowJamNoiOutVec(k) = mean(abs(yOutallMat(:,k)-SigOutallMat(:,k)).^2);
end;
SNROutVec    = 10*log10(PowSigOutVec./PowNoiOutVec);%输出信噪比
SINROutVec   = 10*log10(PowSigOutVec./PowJamNoiOutVec);%输出信干噪比
SINRImpVec   = SINROutVec-SINRInVec;%输出信干噪比改善
JSROutVec    = 10*log10(PowJamOutVec./PowSigOutVec);%输出干信比
JSRImpVec    = JSROutVec-JSRInVec;%输出信干噪比改善
%%%%-------------------------------------------------------------------%%%%

%%    4.4 对有用信号进行码捕获 ------------%%%
disp('    ——对有用信号进行码捕获');
PRNCapFuncMat  = zeros([Len_PRNSnap,1]);  
SigCapFuncMat  = zeros([Len_PRNSnap,1]);    % 用于记录阵列输出信号的码捕获函数
RecCapFuncMat  = zeros([Len_PRNSnap,1]);
yOutCapFuncMat = zeros([Len_PRNSnap,1]);
RelaPeakPosiVec = zeros([Num_Sig,1]);       % 记录各个GPS信号的相关峰位置,也
                                            % 就是其相对延迟
% if FlagDisp==1,
%     figure;     % 后面要在一个界面上画多个图
% end;
for k=1:Num_Sig,
    PRNSnap_k_Vec     = kron(PRNMat(:,k),ones([R,1]));
                                            % 对GPS伪随机码的采样快拍矢量。
    for n=1:Len_PRNSnap,
        PRN_k_ShiftVec       = circshift(PRNSnap_k_Vec,SVTaoRatioVec(k)-(n-1)); %SVTaoRatioVec各个卫星相对于理想零点的延迟
        SigOut_k_Vec_partMat = SigOutallMat((n-1+1):(n-1+Len_PRNSnap),k);   %SigOutallMat:STAP滤波以后的目标信号
        RecSnapVec_partMat   =   RecSnapMat((n-1+1):(n-1+Len_PRNSnap),1);   %RecSnapMat：STAP滤波以前的接收信号
        yOutVec_partMat      =   yOutallMat((n-1+1):(n-1+Len_PRNSnap),k);   %yOutallMat：STAP滤波以后的接收信号
        PRNCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* PRN_k_ShiftVec      /Len_PRNSnap)^2;
        SigCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* SigOut_k_Vec_partMat/Len_PRNSnap).^2;
        RecCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* RecSnapVec_partMat  /Len_PRNSnap).^2;
        yOutCapFuncMat(n,k)  = abs(PRNSnap_k_Vec'* yOutVec_partMat     /Len_PRNSnap).^2;
    end;
    clear SigOut_k_Vec_partMat RecSnapVec_partMat yOutVec_partMat
    PRNCapFuncMat(:,k)  = PRNCapFuncMat(:,k) /max(PRNCapFuncMat(:,k));%归一化？
    SigCapFuncMat(:,k)  = SigCapFuncMat(:,k) /max(SigCapFuncMat(:,k));
    RecCapFuncMat(:,k)  = RecCapFuncMat(:,k) /max(RecCapFuncMat(:,k));
    yOutCapFuncMat(:,k) = yOutCapFuncMat(:,k)/max(yOutCapFuncMat(:,k));
    RelaPeak = max(SigCapFuncMat(:,k));
    RelaPeakPosiVec(k) = find(SigCapFuncMat(:,k)==RelaPeak)-1;% 找到最大相关峰的位置
end;
% if FlagDisp==1,
%     figure;     % 后面要在一个界面上画多个图
%     
%     for k=1:Num_Sig,
%         subplot(2,2,k);
%         tempVec15 = (max(SVTaoRatioVec(k)+1-2*Ppie,1):1:min(SVTaoRatioVec(k)+1+2*Ppie,Len_PRNSnap)).';
%         plot(tempVec15-1,RecCapFuncMat(tempVec15,k), 'b.-');
%         hold on
%         plot(tempVec15-1,PRNCapFuncMat(tempVec15,k), 'k.-');
%         plot(tempVec15-1,yOutCapFuncMat(tempVec15,k), 'g.-');
%         hold off
%         axis([min(tempVec15)-1 max(tempVec15)-1 0 1]);
%         title(sprintf('针对第 %d 个卫星信号的码捕获函数',k));
%     end;   
%    ylabel('相关峰值'); xlabel('码片偏移量(绿:STAP后,蓝:STAP前,黑:理想PRN)');
% end;
%%    4.5 计算权值矢量的方向图 ---------------------------------%%%
    disp(sprintf(' 5、开始计算权值矢量的方向图并画图......'))
    WeightMat = reshape(wVec,M,Ppie);%权值
    switch TypeArray
        case { 'ULA' }, % 针对均匀线阵的仿真结果
            FreqVec = linspace(-0.5*Fsam/Q,0.5*Fsam/Q,201).';
            SitaVec = linspace(-0.5*pi,0.5*pi,201).';
            GainMat = zeros([length(FreqVec),length(SitaVec)]);
            for m=1:M,
                for p=1:Ppie,
                    CosVec = sin(SitaVec);
                    FeqRatioVec = 1+FreqVec/F0;
                    GainSpaceMat = exp(-j*pi*(m-1)*(FeqRatioVec*(CosVec.')));
                    GainTimeMat = kron( exp(-j*2*pi*(p-1)*Q*(F0/Fsam)*FeqRatioVec), ...
                                     ones(1,length(SitaVec)) );
                    GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                end;
            end;
        case { 'NUCA','UCA','UCA_1','UCA_2' }  % 针对加芯均匀圆阵的仿真结果
            % 为了便于演示，固定干扰和信号的Gama均相同，并且也只针对该俯仰角计算空频响应二维图。
            Gama = 45*pi/180; %俯仰角固定
            FreqVec = linspace(-0.5*Fsam/Q,0.5*Fsam/Q,401).';
            SitaVec = linspace(-pi,pi,401).';
            GainMat = zeros([length(FreqVec),length(SitaVec)]);
            switch TypeArray
                case {'NUCA'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-2)/(M-1)-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % 射频频率相对于GPS中心频率的比例
                            GainSpaceMat = (m==1)*ones(size(GainMat)) + ...
                                           (m~=1)*exp(-j*pi*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % 由于阵元位置引入的相差对应的频域响应
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % 由于通道延迟引入的相差对应的频域响应
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case {'UCA'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-1)/M-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % 射频频率相对于GPS中心频率的比例
                            GainSpaceMat = exp(-j*pi*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % 由于阵元位置引入的相差对应的频域响应
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % 由于通道延迟引入的相差对应的频域响应
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case {'UCA_2'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-1)/M-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % 射频频率相对于GPS中心频率的比例
                            Alpha_UCA_2 = 1/sqrt(2);    % 假定相邻阵元的间距为半波长
                            GainSpaceMat = exp(-j*pi*Alpha_UCA_2*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % 由于阵元位置引入的相差对应的频域响应
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % 由于通道延迟引入的相差对应的频域响应
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case { 'UCA_1' }  % 针对均匀圆阵，且以阵元1为参考零点的仿真结果
                    disp('/// ！！！！！！尚未开始编程');
                    GainMat = ones([length(FreqVec),length(SitaVec)]);
            end;%switch2的结束
    end;%switch1的结束
    GainMat_abs = abs(GainMat);
    GainMat_dB = 20*log10(GainMat_abs);
    GainMat_dB_Norm = GainMat_dB - max(max(GainMat_dB));

    if Ppie~=1, % 若时域只有单抽头，则仅画出方向图，否则画出方向频率图
        figure;
        meshc(SitaVec*180/pi,FreqVec/Rc,GainMat_dB_Norm);
        xlabel('角度');  ylabel('归一化频率');
        figure;
        contour(SitaVec*180/pi,FreqVec/Rc,GainMat_dB_Norm,10);
    else
        Freq0Posi = find(FreqVec == 0);
        GainVec = GainMat(Freq0Posi,:).';
        GainVec_dB = 20*log10(abs(GainVec));
        GainVec_dB_Norm = GainVec_dB - max(GainVec_dB);
        figure;
        plot(SitaVec*180/pi,GainVec_dB,'b.-');
        title('P=1,画出方向增益图。');
        hold on;
        for i=1:Num_Jam
            yer=linspace(5,-70,201);
            xer=ones(201,1)*Jam(i);
            plot(xer,yer);
        end
%         tempMat = kron(GainVec_dB_Norm,ones([1,length(FreqVec)]));
        tempMat = kron(ones([length(FreqVec),1]),GainVec_dB.');
%         figure;
%         meshc(SitaVec*180/pi,FreqVec/Rc,tempMat);
%         xlabel('角度');  ylabel('归一化频率');
%         figure;
%         contour(SitaVec*180/pi,FreqVec/Rc,tempMat,10);
    end;
%% % 程序结束
Time_end = cputime;
Time_used = Time_end - Time_begin;    % 整个仿真过程耗费的时间
Time_clock = fix(clock);
disp(sprintf('程序结束时间： %d-%d-%d  %d:%d:%d',...
              Time_clock(1),Time_clock(2),Time_clock(3), ...
              Time_clock(4),Time_clock(5),Time_clock(6)));
disp(sprintf('整个仿真程序所用时间： Time_used = %g 秒。',Time_used));
disp(sprintf('------程序结束------'));
disp(sprintf('  '));
disp(sprintf('  '));


%1.功率倒置算法很简单，已知导向矢量，直接用MVDR准则就行
%2.经典维纳滤波是经过导向矢量求出期望信号，进而求出阵列接收数据与期望信号的互相关矢量
%   再用维纳解求出最优权矢量。
%3.其他多级维纳滤波算法则是在经典维纳滤波基础上就行降秩处理。

