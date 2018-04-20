% function SINROutVec = ClassicalWeinerForKongpin(SigSitaVec,SigGamaVec,P)
%master bench   %696
clear all;
clc;
format long e
% ����ʼ����ʾ����ʼʱ�䡣
disp(sprintf('------����ʼ���� ......'));
%%%%%%%%%%%%%%%%%%%%%%%%%%

SigSitaVec = (pi/180)*[10];       % �ź�����
Jam=[50;10;-30];
JamSitaVec = Jam*pi/180; %���ŷ�λ��  ƽ������
AlgorithmType = '5';                                      %ѡ���㷨ִ������
ybsn=0.0001;
C=1;
dirta=2;
P=120;

Num_Sig =size(SigSitaVec);                                                                 %GPS�����ź���Ŀ
Num_Jam = size(Jam);                                                              %���Ÿ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_begin = cputime;
Time_clock = fix(clock);
disp(sprintf('����ʼʱ�䣺 %d-%d-%d  %d:%d:%d',...
              Time_clock(1),Time_clock(2),Time_clock(3), ...
              Time_clock(4),Time_clock(5),Time_clock(6)));
FlagDisp = 1;
%FlagFixLen = 1;
TypeArray = 'ULA';  %%%����������״ѡ����Բ�ͬ�������ɵ���ʸ�������SteerVecGen��������
%  SttrandBegin =round(rand*100000);
SttrandBegin = 70;
rand('state',SttrandBegin);
SttrandConfg = round(rand*10000);%����һ��4λ���������
SttrandBit   = round(rand*10000);
SttrandnJam  = round(rand*10000);   %���ţ�    
SttrandnNoi  = round(rand*10000);   %������

%JamSitaVec = [70;90;150;130]*pi/180; %���ŷ�λ��   ��������
% JamSitaVec = [-30;30;70]*pi/180; %���ŷ�λ��  ƽ������
JamGamaVec = [45;45; 45]*pi/180;%���Ÿ�����
% SigSitaVec = 10;%�źŷ�λ
% SigGamaVec = 30;%�źŸ���

disp(sprintf(' 1�����������趨 '));
%% % 1�����������趨 --------------------------------------------%%%%%%
%% %    1.1 �������к����Ͽ�ʱ�����ʱ���ͷ���� --------------------%%%%
rand('state',SttrandConfg);            % װ������������α�������������״̬,�磺ִ��rand('state',0);rand(1)   rand('state',SttrandConfg);rand(1)
M = 10;                  % ��Ԫ����
Ppie =1;                  %ʱ���ͷ����
if ((Ppie-1)/2-round((Ppie-1)/2))>10^(-12),
    disp('���澯��Ҫ��ʱ���ͷ���� P ����������');
else
    P_h = round((Ppie+1)/2);
    DlySTAP = P_h - 1;  % STAP�������ӳ�,��Ȼ���Դ���[0��P-1]֮�䣬��Ŀǰ�̶�Ϊ2��%�ⲻ̫��ΪʲôҪ��ô��
end;

%%    1.2 ����Ƶ�ʡ�������������ȡ���� ---------------------------%%%
Rc = 1.023*10^6;        % GPS C/A�����Ƭ����
F0 = 1540*Rc;           % GPS�źŵ�L1Ƶ�ʣ���ƵƵ��
Fc =   70*Rc;           % �ز�Ƶ�ʣ����ô�ͨ��������
G  =   40*2;              % ÿһ����Ƭ�����ڵĲ���λ����Ҫ�����Ƶ�ʱ�������Ƭ���ʵ�ż������
G_h = G/2;
G_q = G/4;
Q  =    8;              % ��ȡ����
R  = G / Q;             % ��ȡ��ÿ����Ƭʱ�������ڵĲ�����������������
R_h = R/2;
Fsam =  G*Rc;           % Ҫ�����Ƶ������Ƭ���ʵ�ż������
Num_Bit = 2;            % ��Ϣ���ص���Ŀ
Len_PRN = 1023;         % ÿ����Ƶ�볤1023����Ƭ
Len_PRNSnap = Len_PRN*R;% ÿ��������PRN�����еĲ���������
Len_PRNSam  = Len_PRN*G;% ÿ��������PRN�����еĲ�������
Num_Chip = Num_Bit  * Len_PRN;      % ��Ƭ����
Num_Sam  = Num_Chip * G;            % �������еĳ���
Num_Snap = Num_Chip * R;            % ���Ͽ�ʱ�����Ŀ��������������д������
                                    % �����ݳ���
%%    1.3 ���ջ����ߴ�����������(����) ---------------------------%%%
NoiPow = 1;                         % �������Ĺ���
%%    1.4 GPS�źŵĹ��ʡ����Ⱥͳ��� -----------------------------%%%

switch TypeArray                    % ������������źŵķ�λ�Ǻ͸�����
    case {'ULA'}
%         SigSitaVec = (pi/180)*[0];% �������ǵķ�λ��
        SigGamaVec = (pi/180)*[45];% �������ǵĸ����ǣ���������
                                                        % ��ȡֵ������
    case {'NUCA','UCA','UCA_1','UCA_2'}
        SigSitaVec = (pi/180)*[ 150; 170;30;200;251;307];% �������ǵķ�λ��
        SigGamaVec = (pi/180)*[ 45; 45; 45; 45; 45; 45];% �������ǵĸ�����
end;
SVIDVec = [1;5;9;13;17;19];               % �����������
SVTaoRatioVec = [20+round(50*rand([Num_Sig,1]))];   % ����������������������ӳ�/Tsnap��    round()��������
SNRVec = -27.5 * ones([Num_Sig,1]); % ���������źŵ������(dBֵ)��ȡΪGPS   ones([Num_Sig,1])�źŸ����У��м����źž��м��������
                                    % C/A��ĵ�������ȣ�������7dB��Ϊ���ʺ�
                                    % �����STAPǰSNR��ֵ�ﵽ-21.9dB��                                    
SigPowVec = NoiPow*10;% ����GPS�źŵĽ��չ��ʡ�������Ⱥ��������ʣ������źŹ���

SigFaiVec = rand([Num_Sig,1])*2*pi; % �����ź��ز��ĳ��࣬���ڽ��ջ�����������
                                    % ���ĳ����ǹ̶��ģ���˿�����Ϊ�����źŵ�    ��̫���ף���������
                                    % �ز�����������ֲ��ġ�
                                   
%%    1.5 ���ŵ���Ŀ�����͡����ʺ������ -------------------------%%%
 JamKindVec = ['FBJam';'FBJam';'FBJam'];    % ��������
%JamKindVec = ['FBJam';'FBJam';'FBJam';'STJam'];    % ��������
                                       % ='FBJam' ����������ţ�
                                       % ='NBJam' խ���������ţ�
                                       % ='STJam' ��Ƶ���ţ�
JNRVec = [30;30;30];                % ����ȣ���������������Ĺ��ʱ�(dBֵ)
%JNRVec = [40;40;35;35;35]-7.5;        % ����ȣ���������������Ĺ��ʱ�(dBֵ)
%JamPowVec = NoiPow*10.^(JNRVec/10);    % �������ŷ����Ĺ���(����ָ���Ǹ��ŵĸ���
JamPowVec = NoiPow*10.^(JNRVec/10);          % ��,ʵ�����Ǹ����ʵ�һ��)  
                                       
JamFrqOffVec = [0; -1; 0.5;-1;+1]*Rc;          
                                       % ��������Ƶ��������ز���Ƶ�ʲ�,������Ƶ
                                       % ���������塣�䶨��Ϊ(fi-Fc)��FcΪ�ź���
                                       % �������ڻ������������Fc=0��
                                       % ���ǵ���Ƶ��ͨ�˲�����ͨ����Χ����˸���
                                       % Ƶ��ϵ���ľ���ֵ���ܳ���1��
JamFaiVec = rand([Num_Jam,1])*2*pi;    % ���ŵ��ز�����

%%    1.6 ��ͨ�˲���,����һ��ͨ������Ƶ����0.25*Fsam�Ĵ�ͨ�˲��� ---%%%
Rp = 1;                             % ͨ����������Ʋ�
Rs = 80;                            % ���˥��
% f1 = [0.2/15.6  6.5/15.6  9.1/15.6  15.4/15.6];    % ͨ�������ֹƵ���趨
f1 = [0.1/G_h  (G_q-1)/G_h  (G_q+1)/G_h  (G_h-0.1)/G_h];
                                    % ͨ�������ֹƵ���趨
a1 = [    0    1    0     ];        % ָ��ͨ�������
dev1 = [10^(-Rs/20)  (10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)]; 
                                    % ����ͨ���Ʋ������˥��׼������
[n,fo,ao,wo] = firpmord(f1,a1,dev1); % ȷ��FIR�˲����Ľ���
n1 = ceil(n/2)*2 + 10;      % ����firpm.m���FIR�˲�����������΢�ﲻ�����
                            % ˥��Ҫ�󣬹���ǰ���һЩ���ȡ����ң�Ҫ��n1��
                            % ż�����Ա�����Ⱥ�ӳ�n1/2��һ��������
BPFVec = firpm(n1,fo,ao,wo); % ����������λ�Ĵ�ͨ�˲���
GrpDlyBPF = round((n1-1)/2);    % ���Դ�ͨ�˲�����Ⱥ�ӳ�
%%%%--------------------------------------------------------------------%%%
%%    1.7 ��ͨ�˲������ ---------------------------------------%%%
%            ���������±�ƵʱI��Q��·��Ƶ�����ݵĵ�ͨ�˲�����׼����ȡ
Rp = 1;                             % ͨ����������Ʋ�
Rs = 50;                            % ���˥��
% f2 = [1.1/15.6  2.5/15.6];           % ͨ�������ֹƵ���趨
f2 = [1/G_h  R_h/G_h];                % ͨ�������ֹƵ���趨
a2 = [ 1    0  ];                   % ָ��ͨ�������
dev2 = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)]; 
                                    % ����ͨ���Ʋ������˥��׼������
[n,fo,ao,wo] = firpmord(f2,a2,dev2); % ȷ��FIR�˲����Ľ���
n2 = ceil(n/2)*2 + 4;%10;      % ����firpm.m���FIR�˲�����������΢�ﲻ�����
                            % ˥��Ҫ�󣬹���ǰ���һЩ���ȡ����ң�Ҫ��n2��
                            % ż�����Ա�����Ⱥ�ӳ�n2/2��һ��������
LPFVec = firpm(n2,fo,ao,wo);         % ����������λ�ĵ�ͨ�˲���
GrpDlyLPF = round((n2-1)/2);        % ���Ե�ͨ�˲�����Ⱥ�ӳ�

FlagLoJit = 0;                      % ����ʱ����ʱ�򶶶�����
% LoJit_std = 2*10^(-12);             % ����Ƶ�ʺ�ADC����Ƶ�ʵĺϳɵ�Ч
%                                     % ʱ�򶶶��ķ���
% if FlagLoJit == 1,
%     disp(sprintf('  (9)�����п����˱���Ƶ�ʺ�ADC����Ƶ�ʵ�ʱ�򶶶������ߺϳɵ�Ч��ʱ�򶶶��ķ���Ϊ %4.2g ps��',...
%                    LoJit_std*10^12) );
% end;
%%%%-------------------------------------------------------------------%%%%
%%   2����������׼�� ------------------------------------------%%%%
 disp(sprintf(' 2������Ƶ��ʼ�������ݣ�����Ƶ�͵�ͨ�˲�����ȡ��õ��������ݡ�'));
%%     2.1 �����������˵�� ------------------------------------%%%
        Num_BitRaw = Num_Bit + ceil((GrpDlyBPF+GrpDlyLPF)/Len_PRNSam/2)*2;
                                                % ԭʼ��Ϣ���ظ�����Ҫ����ż��
        Num_ChipRaw = Num_BitRaw * Len_PRN;     % ԭʼ��Ƭ����
        Num_SamRaw  = Num_ChipRaw * G;          % ԭʼ��������
        Num_SamADC = Num_SamRaw - GrpDlyBPF;    % ADC�����õ����������г���
        nVec = ( 1:1:Num_SamRaw ).';            % ���������ʸ������1��ʼ 
%%     2.2 ���C/A�����������Ƭ���ݺ���Ƭ�������� ----------------%%%
        Ratio_Sig  = Fc/Fsam;                   %һ����Ƭ�����㣬��Ҫ�����ز�����
        BitMat_Raw = zeros([Num_BitRaw,Num_Sig]);
        PRNMat     = zeros([ Len_PRN  ,Num_Sig]);
        DSSnapMat  = zeros([round(Num_SamRaw/Q) ,Num_Sig]);
        SigVecMat  = zeros([Num_SamADC,Num_Sig]);
        SigIFallMat= zeros([Num_SamADC,   M   ,Num_Sig]);%M = 4; % ��Ԫ���� Num_Sig �źŸ���
        SigIFMat   = zeros([Num_SamADC,   M   ]);
        SigSteerallMat= zeros([M,Num_Sig]);
        for k=1:Num_Sig,
            rand('state',SttrandBit+(k-1)*237); % ���ȷֲ�α�����״̬װ��
            BitMat_Raw(:,k) = [ AJ_AryBalance(rand([           Num_Bit,1]));  ...
                                AJ_AryBalance(rand([Num_BitRaw-Num_Bit,1])) ] ;
                                    % ����һ��+1��-1������ȵľ������У�����
                                    % Ҫ��ǰNBit������Ԫ��֮��Ҳ����0��
            [G1Vec,PRNMat(:,k)] = AJ_GPSCAPRNGen(SVIDVec(k)); 
                                    % ������k��GPS���ǵ�PRN�룬PRNMatΪ����1023��Ƭ��CA��
            DSChipVec = kron( BitMat_Raw(:,k),PRNMat(:,k) );
            DSSnapVec_NoTao = kron(DSChipVec,ones([R,1]));
                                    % ��GPS�ź���Ƭ�Ŀ��Ĳ�������ʸ������δ
                                    % �����ӳ�
            DSSnapMat(:,k) = circshift(DSSnapVec_NoTao,SVTaoRatioVec(k));
                                    % �Ѹ�������PRN�������ӳ����ֳ���
            DSSamVec = kron(DSSnapMat(:,k),ones([Q,1]));
            tempVec1 = exp( j*(2*pi*Ratio_Sig*(nVec-1)+SigFaiVec(k) ...
                                     +0.5*(1- DSSamVec)*pi)        ) ;
                                    % ��Ƭ+1��Ӧ�ڵ�����λ0����Ƭ-1��Ӧ�ڵ���
                                    % ��λpi�������Ժ�ָ����
                                    %tempVec1Ϊģ���±�Ƶ�����Ƶ�źţ�%����SigFaiVec�źų��࣬��ʱ�Ѿ��������Ƶ����
            SigVecMat(:,k) = IFFltandPowSet(tempVec1,SigPowVec(k),BPFVec);
                                    %��ģ���±�Ƶ�����Ƶ�źŽ���ģ���ͨ�˲�������SigPowVec�źŹ���
            SigSteerVec = SteerVecGen(TypeArray,SigSitaVec(k),SigGamaVec(k),M);
            Alb=SigSteerVec;
                                    %��Բ�ͬ�������ɵ���ʸ����TypeArray��Ԫ���ͣ��źŵķ�λ��������Ԫ����
            SigSteerallMat(:,k) = SigSteerVec;
                                    %�����źŸ���Num_Sig����������
            SigIFallMat(:,:,k) = kron(SigVecMat(:,k),SigSteerVec.');
                                    %���е���ʸ������Ƶ�ź����
            SigIFMat = SigIFMat + SigIFallMat(:,:,k);
                                    % �õ�ÿ�����߶������źŵĸ�����������,��
                                    % �����ۼ���������Ϊ�����ź�
        end;
        clear tempVec1
        clear SigVecMat
        
%%     2.3 �������� -------------------------------------------%%%
        randn('state',SttrandnJam);             % ��˹�ֲ�α�����״̬װ��
        JamIFMat = zeros([Num_SamADC,M]);       % ��������յ��ĸ��Ÿ�������
        JamSteerallMat = zeros([M,Num_Jam]);
        for k=1:Num_Jam,
            switch JamKindVec(k,:),
                case {'FBJam'},                 % �����Ǹ�˹����������
                    tempVec1 =    randn(Num_SamRaw+GrpDlyLPF,1) ...
                               +j*randn(Num_SamRaw+GrpDlyLPF,1) ;
                                                % ����һ����������������
                    tempVec0 = filter(LPFVec,1,tempVec1);
                    Ratio_Jam = Ratio_Sig;      % �����������ŵ�����Ƶ�����ź���ͬ��
                    tempVec2 =    tempVec0((GrpDlyLPF+1):(GrpDlyLPF+Num_SamRaw)) ... 
                              .* exp(j*(2*pi*Ratio_Jam*(nVec-1)+JamFaiVec(k)));
                case {'PBJam'},
                  
                case {'STJam'},                 % �����ǵ�Ƶ����
                    Ratio_Jam = (Fc+JamFrqOffVec(k))/Fsam;
                    tempVec2 = exp(j*(2*pi*Ratio_Jam*(nVec-1)+JamFaiVec(k)));
                                                % �������ḳ��
            end;
            JamVec = IFFltandPowSet(tempVec2,JamPowVec(k),BPFVec);
            SteerVec_Jam = SteerVecGen(TypeArray,JamSitaVec(k),JamGamaVec(k),M);
                                            % �������߽��յ����źŸ�����λ��
                                            % ����Ϊ�Խ����źŵ�ģ1������ˡ�
%             JamSteerallMat(:,k) = SteerVec_Jam;
            JamIFMat = JamIFMat + kron(JamVec,SteerVec_Jam.');
                                            % �����п�����ŵĸ����������ݵ�������        
        end;
        clear tempVec2
%%     2.4 ���������� -----------------------------------------%%%
        randn('state',SttrandnNoi);         % ��˹�ֲ�α�����״̬װ��
        NoiIFMat = sqrt(NoiPow/2)*(      randn(Num_SamADC,M) ...
                                     + j*randn(Num_SamADC,M) );
%%    2.5 �ź��������˵�� ------------------------------------%%%
        
%%    2.6 �źš����źͰ������ϳ�ΪADC���������  -----------------%%%
        RecIFMat = SigIFMat + JamIFMat + NoiIFMat;        
                                            % �õ�������ĸ����������ݾ���   
%         if FlagLoJit == 1,                  % �жϷ���ʱ�Ƿ���ʱ�򶶶�����
%             LoJitMat = exp(j*2*pi*F0*LoJit_std*randn(size(RecIFMat)));
%             RecIFMat = RecIFMat.*LoJitMat;
%         end;
      
        RecADCMat = real(RecIFMat);         % ʵ����ADC�����õ�����ʵ������        
        SigADCMat = real(SigIFMat);         %Ŀ���ź�ʵ������
        JamADCMat = real(JamIFMat);         %����ʵ������
        NoiADCMat = real(NoiIFMat);         %����ʵ������
        SigADCallMat = real(SigIFallMat);   %�źŸ���û�е��ӵ�����
        
%%    2.7 ��Ƶ�͵�ͨ�˲�(��Ҫ���ճ�ȡ��Ҫ��ȷ�������ֹƵ��)-------%%%
        Fai_Mix = 0;                    % ���ػ�Ƶ�ز�����
        RecBFMat = MixandLPF(RecADCMat,Ratio_Sig,Fai_Mix,LPFVec);%% ��ADC�����õ��Ĵ�ͨ�������н��л�Ƶ�͵�ͨ�˲������ҽ�I��Q��·����ź��������Ϊһ�������Ļ����ź����С�
        SigBFMat = MixandLPF(SigADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        JamBFMat = MixandLPF(JamADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        NoiBFMat = MixandLPF(NoiADCMat,Ratio_Sig,Fai_Mix,LPFVec);
        SigBFallMat = zeros([Num_SamADC-GrpDlyLPF,M,Num_Sig]);
        for k=1:Num_Sig,
            SigBFallMat(:,:,k) = MixandLPF(SigADCallMat(:,:,k),Ratio_Sig,Fai_Mix,LPFVec);
        end;
%%     2.8 ��ȡ -----------------------------------------------%%%
        nVec3 = ( (0:1:Num_Snap-1)*Q+1 ).'; % Q����ȡ��λ�úŹ���һ��ʸ��
        RecSnapMat = RecBFMat(nVec3,:);     % �����źŻ�Ƶ�˲����Q����ȡ
        SigSnapMat = SigBFMat(nVec3,:);     % �����źŻ�Ƶ�˲����Q����ȡ
        JamSnapMat = JamBFMat(nVec3,:);     % ���Ż�Ƶ�˲����Q����ȡ
        NoiSnapMat = NoiBFMat(nVec3,:);     % ������Ƶ�˲����Q����ȡ
        SigSnapallMat = SigBFallMat(nVec3,:,:);
%% %----------------------�����ź��������-----------------------%%%
%%  3������Ӧ�˲� ---------------------------------------------%%%%
%%    3.1 ���������趨������׼�� --------------------------------%%%
disp(' 3����ʼ���Ȩֵʸ��');
RecSTAPMat = zeros([M*Ppie,Num_Snap-(Ppie-1)]);
% ��x(m,n)��ʾ��m������֧·�ϵ�n��ʱ�̵Ŀ���ֵ�����n��ʱ��STAP�˲���ʸ���涨Ϊ��
% [x(1,n+P-1); ... x(M,n+P-1);x(1,n+P-2); ... x(M,n+P-2); ...;x(1,n); ...
% x(M,n)]       
for n=1:(Num_Snap-(Ppie-1)),
    RecSTAPMat(:,n) = reshape(RecSnapMat((n+Ppie-1):-1:n,:).',M*Ppie,1); %��ʵ�Ľ����ź�����
    SigSTAPMat(:,n) = reshape(SigSnapMat((n+Ppie-1):-1:n,:).',M*Ppie,1); %��һ����Ҫ�õ�����Ŀ���ź�����   
end;

%%    3.2 �˲�Ȩֵ���� ----------------------------------------%%%
% AlgorithmType = '5';                                      %ѡ���㷨ִ������
%'0'�����ʵ����㷨
%'1'��MVDR�㷨
%'2'������ά���˲��㷨
%'3'��CSA_MSNWF:�������ṹ�Ķ༶ά���˲�    P72
%'4'���ಽ���������㷨    P73
disp(sprintf('----ִ�е�%s���㷨 ',AlgorithmType));
Time_begin1 = cputime;
%%%%%%%%%%%%%ר�����ڼ���Rx%%%%%%%%%%%%%%%%%%
        L =1024;      
        xAnaMat = RecSTAPMat(:,1:L);%�����ź�
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %��������������ؾ���
        Rx=RxMat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch AlgorithmType, 
      %% ���ʵ����㷨  
        case {'0'},
        disp('    �������ù��ʵ����㷨����Ȩֵ����,');
        L =1024;
        xAnaMat = RecSTAPMat(:,1:L);%�����ź�
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %��������������ؾ���
        Rx=RxMat;
        a=zeros([M*Ppie,1]);
        a(1)=1;  % ���ʵ���
        wVec = (inv(RxMat)*a)/(a'*inv(RxMat)*a);
  %%  MVDR�㷨
        case {'1'},
        disp('    ��������MVDR�㷨����Ȩֵ����,');
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);%�����ź�
        [MM,NumWeiner] = size(xAnaMat);
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %��������������ؾ���
        a=zeros([M*Ppie,Num_Sig]);
        %a=st_cir(SigSteerVec,P,Fsam,Fc); %��һ���ǵ���ʸ��
        for k=1:Num_Sig
             a(:,k)=st_cir(SigSteerallMat(:,k),Ppie,Fsam,Fc);
        end     
        I_sig=ones(Num_Sig,1);
        wVec = (inv(RxMat)*a)/(a'*inv(RxMat)*a)*I_sig;
        %% %-------------------��1�־���ά���˲�-------------------------%%
    case {'2'},
        disp('    �������þ���ά���˲��㷨����Ȩֵ����,');
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);%�����ź�
        [MM,NumWeiner] = size(xAnaMat);
        dRow    = SigSTAPMat( 1 ,1:L);%�����ź�
      
        RxMat = xAnaMat*xAnaMat'/NumWeiner;     %��������������ؾ���
        pxVec = xAnaMat*dRow'   /NumWeiner;  %�����ź��������źŵĻ����ʸ��
        wVec = zeros([MM,1]);
        wVec =   (RxMat \ pxVec); % ���þ��������ͨ����˹��Ԫ�����ά�����Ž⡣���ʽA\B�����൱��inv��A��*B
        wVec = wVec/norm(wVec);
%% %----------------------%��2��CSA_MSNWF_wVec:�������ṹ�Ķ༶ά���˲���%----%
    case {'3'},
%         TypeConst = 'SimpConst';
        disp('    ��������CSA_MSWNF�㷨����Ȩֵ����,');
        L = 1024;    % ��L����XX'��ƽ��������hi��wi�����
                                        % ������128�����ݵ㡣Ҫ��L��2������
                                        % ���ݣ����������ƽ��ֵ�Ĳ��������
                                        % ��
        xAnaMat = RecSTAPMat(:,1:L);
        TypeConst = 'SSSCConst';
        wVec = CSA_MSNWF_wVec(M,Ppie,L,TypeConst,SigSteerVec,xAnaMat);                               
                                        

%% %-----------------------------��3�֣�ICSA_MSNWF_wVec:�ಽ���������������㷨--------------------------%%
    case {'4'},
        L = 1024;
        xAnaMat = RecSTAPMat(:,1:L);
        TypeConst = 'SSSCConst';
        wVec = ICSA_MSNWF_wVec(M,Ppie,L,TypeConst,SigSteerallMat(:,1),xAnaMat);
    case{'5'}
        Kmax=9999; %maxminum loop time
        minchaju=0.1;%the minum varience to judge if it is regressed. warining!
        bi=1;   %warning! used in caculation J2 method in function jhs 
        g=1; % ���������ķ�Ӧ
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
% ����i=1;
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
end      %�㷨����ѡ�������


        wVecallMat = kron(ones([1,Num_Sig]),wVec);  
%% %--------------------------------------------------------%%
Time_mid = cputime;
Time_used = Time_mid - Time_begin1;    
disp(sprintf('----��Ȩ����ʱ�䣺 Time_used = %g �롣',Time_used));
%% %%%% 4���˲�,�����������ź� ------------------------------------------%%%%%%
%% %    4.1 �˲� -----------------------------------------------%%%%
disp(' 4��Ȩֵʸ��������ϣ����ڿ�ʼ�ø�Ȩֵʸ�����п������˲�');
yOutallMat   = zeros([Num_Snap,Num_Sig]);%NumSnap��ȡ������ݵ���
JamOutallMat = zeros([Num_Snap,Num_Sig]);%����
NoiOutallMat = zeros([Num_Snap,Num_Sig]);%����
SigOutallMat = zeros([Num_Snap,Num_Sig]);%�����ź�
for k=1:Num_Sig
    yOutallMat(:,k)   = STAPFilter(RecSnapMat,Ppie,wVec,DlySTAP);%��ʵ���ݣ�DlySTAP��stap�ӳ�
    SigOutallMat(:,k) = STAPFilter(SigSnapallMat(:,:,k),Ppie,wVec,DlySTAP);%�����ź�
    JamOutallMat(:,k) = STAPFilter(JamSnapMat,          Ppie,wVec,DlySTAP);%����
    NoiOutallMat(:,k) = STAPFilter(NoiSnapMat,          Ppie,wVec,DlySTAP);%����
end;
%% %--------------------------------------------------------%%
Time_mid = cputime;
Time_used = Time_mid - Time_begin1;    
disp(sprintf('----�˲�����ʱ�䣺 Time_used = %g �롣',Time_used));
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %    4.2 ����STAP�����STAP���SINR ---------------------------%%%%
%%% ������������£����ÿһ�����ߵ����Ȩֵ�����õ�����������Ÿ����
%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowSigInVec    = zeros([Num_Sig,1]);
PowJamInVec    = zeros([Num_Sig,1]);
PowNoiInVec    = zeros([Num_Sig,1]);
PowJamNoiInVec = zeros([Num_Sig,1]);
for k=1:Num_Sig,
    PowSigInVec(k)    = mean(abs(               SigSnapallMat(:,1,k)).^2);%�����źŹ���
    PowJamInVec(k)    = mean(abs(               JamSnapMat(:,1)    ).^2);%������Ź���
    PowNoiInVec(k)    = mean(abs(               NoiSnapMat(:,1)    ).^2);%������������
    PowJamNoiInVec(k) = mean(abs(RecSnapMat(:,1)-SigSnapallMat(:,1,k)).^2);%������ż���������

end;
SNRInVec     = 10*log10(PowSigInVec./PowNoiInVec);%���������
SINRInVec    = 10*log10(PowSigInVec./PowJamNoiInVec);%%�����Ÿ����
JSRInVec     = 10*log10(PowJamInVec./PowSigInVec);%������ű�
%%%%%%%%%%%%%%%%���%%%%%%%%%%%%%%%%%%%%%%%%
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
SNROutVec    = 10*log10(PowSigOutVec./PowNoiOutVec);%��������
SINROutVec   = 10*log10(PowSigOutVec./PowJamNoiOutVec);%����Ÿ����
SINRImpVec   = SINROutVec-SINRInVec;%����Ÿ���ȸ���
JSROutVec    = 10*log10(PowJamOutVec./PowSigOutVec);%������ű�
JSRImpVec    = JSROutVec-JSRInVec;%����Ÿ���ȸ���
%%%%-------------------------------------------------------------------%%%%

%%    4.4 �������źŽ����벶�� ------------%%%
disp('    �����������źŽ����벶��');
PRNCapFuncMat  = zeros([Len_PRNSnap,1]);  
SigCapFuncMat  = zeros([Len_PRNSnap,1]);    % ���ڼ�¼��������źŵ��벶����
RecCapFuncMat  = zeros([Len_PRNSnap,1]);
yOutCapFuncMat = zeros([Len_PRNSnap,1]);
RelaPeakPosiVec = zeros([Num_Sig,1]);       % ��¼����GPS�źŵ���ط�λ��,Ҳ
                                            % ����������ӳ�
% if FlagDisp==1,
%     figure;     % ����Ҫ��һ�������ϻ����ͼ
% end;
for k=1:Num_Sig,
    PRNSnap_k_Vec     = kron(PRNMat(:,k),ones([R,1]));
                                            % ��GPSα�����Ĳ�������ʸ����
    for n=1:Len_PRNSnap,
        PRN_k_ShiftVec       = circshift(PRNSnap_k_Vec,SVTaoRatioVec(k)-(n-1)); %SVTaoRatioVec����������������������ӳ�
        SigOut_k_Vec_partMat = SigOutallMat((n-1+1):(n-1+Len_PRNSnap),k);   %SigOutallMat:STAP�˲��Ժ��Ŀ���ź�
        RecSnapVec_partMat   =   RecSnapMat((n-1+1):(n-1+Len_PRNSnap),1);   %RecSnapMat��STAP�˲���ǰ�Ľ����ź�
        yOutVec_partMat      =   yOutallMat((n-1+1):(n-1+Len_PRNSnap),k);   %yOutallMat��STAP�˲��Ժ�Ľ����ź�
        PRNCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* PRN_k_ShiftVec      /Len_PRNSnap)^2;
        SigCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* SigOut_k_Vec_partMat/Len_PRNSnap).^2;
        RecCapFuncMat(n,k)   = abs(PRNSnap_k_Vec'* RecSnapVec_partMat  /Len_PRNSnap).^2;
        yOutCapFuncMat(n,k)  = abs(PRNSnap_k_Vec'* yOutVec_partMat     /Len_PRNSnap).^2;
    end;
    clear SigOut_k_Vec_partMat RecSnapVec_partMat yOutVec_partMat
    PRNCapFuncMat(:,k)  = PRNCapFuncMat(:,k) /max(PRNCapFuncMat(:,k));%��һ����
    SigCapFuncMat(:,k)  = SigCapFuncMat(:,k) /max(SigCapFuncMat(:,k));
    RecCapFuncMat(:,k)  = RecCapFuncMat(:,k) /max(RecCapFuncMat(:,k));
    yOutCapFuncMat(:,k) = yOutCapFuncMat(:,k)/max(yOutCapFuncMat(:,k));
    RelaPeak = max(SigCapFuncMat(:,k));
    RelaPeakPosiVec(k) = find(SigCapFuncMat(:,k)==RelaPeak)-1;% �ҵ������ط��λ��
end;
% if FlagDisp==1,
%     figure;     % ����Ҫ��һ�������ϻ����ͼ
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
%         title(sprintf('��Ե� %d �������źŵ��벶����',k));
%     end;   
%    ylabel('��ط�ֵ'); xlabel('��Ƭƫ����(��:STAP��,��:STAPǰ,��:����PRN)');
% end;
%%    4.5 ����Ȩֵʸ���ķ���ͼ ---------------------------------%%%
    disp(sprintf(' 5����ʼ����Ȩֵʸ���ķ���ͼ����ͼ......'))
    WeightMat = reshape(wVec,M,Ppie);%Ȩֵ
    switch TypeArray
        case { 'ULA' }, % ��Ծ�������ķ�����
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
        case { 'NUCA','UCA','UCA_1','UCA_2' }  % ��Լ�о����Բ��ķ�����
            % Ϊ�˱�����ʾ���̶����ź��źŵ�Gama����ͬ������Ҳֻ��Ըø����Ǽ����Ƶ��Ӧ��άͼ��
            Gama = 45*pi/180; %�����ǹ̶�
            FreqVec = linspace(-0.5*Fsam/Q,0.5*Fsam/Q,401).';
            SitaVec = linspace(-pi,pi,401).';
            GainMat = zeros([length(FreqVec),length(SitaVec)]);
            switch TypeArray
                case {'NUCA'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-2)/(M-1)-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % ��ƵƵ�������GPS����Ƶ�ʵı���
                            GainSpaceMat = (m==1)*ones(size(GainMat)) + ...
                                           (m~=1)*exp(-j*pi*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % ������Ԫλ�����������Ӧ��Ƶ����Ӧ
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % ����ͨ���ӳ����������Ӧ��Ƶ����Ӧ
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case {'UCA'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-1)/M-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % ��ƵƵ�������GPS����Ƶ�ʵı���
                            GainSpaceMat = exp(-j*pi*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % ������Ԫλ�����������Ӧ��Ƶ����Ӧ
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % ����ͨ���ӳ����������Ӧ��Ƶ����Ӧ
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case {'UCA_2'},
                    for m=1:M,
                        for p=1:Ppie,
                            CosVec = cos(2*pi*(m-1)/M-SitaVec);
                            FeqRatioVec = 1+FreqVec/F0; % ��ƵƵ�������GPS����Ƶ�ʵı���
                            Alpha_UCA_2 = 1/sqrt(2);    % �ٶ�������Ԫ�ļ��Ϊ�벨��
                            GainSpaceMat = exp(-j*pi*Alpha_UCA_2*cos(Gama)*(FeqRatioVec*(CosVec.')));
                                                        % ������Ԫλ�����������Ӧ��Ƶ����Ӧ
                            GainTimeMat = kron( exp(-j*2*pi*(p-1)*F0*(Q/Fsam)*FeqRatioVec), ...
                                                ones(1,length(SitaVec)) );
                                                        % ����ͨ���ӳ����������Ӧ��Ƶ����Ӧ
                            GainMat = GainMat + WeightMat(m,p)'*(GainTimeMat.*GainSpaceMat);
                        end;
                    end;
                case { 'UCA_1' }  % ��Ծ���Բ��������Ԫ1Ϊ�ο����ķ�����
                    disp('/// ��������������δ��ʼ���');
                    GainMat = ones([length(FreqVec),length(SitaVec)]);
            end;%switch2�Ľ���
    end;%switch1�Ľ���
    GainMat_abs = abs(GainMat);
    GainMat_dB = 20*log10(GainMat_abs);
    GainMat_dB_Norm = GainMat_dB - max(max(GainMat_dB));

    if Ppie~=1, % ��ʱ��ֻ�е���ͷ�������������ͼ�����򻭳�����Ƶ��ͼ
        figure;
        meshc(SitaVec*180/pi,FreqVec/Rc,GainMat_dB_Norm);
        xlabel('�Ƕ�');  ylabel('��һ��Ƶ��');
        figure;
        contour(SitaVec*180/pi,FreqVec/Rc,GainMat_dB_Norm,10);
    else
        Freq0Posi = find(FreqVec == 0);
        GainVec = GainMat(Freq0Posi,:).';
        GainVec_dB = 20*log10(abs(GainVec));
        GainVec_dB_Norm = GainVec_dB - max(GainVec_dB);
        figure;
        plot(SitaVec*180/pi,GainVec_dB,'b.-');
        title('P=1,������������ͼ��');
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
%         xlabel('�Ƕ�');  ylabel('��һ��Ƶ��');
%         figure;
%         contour(SitaVec*180/pi,FreqVec/Rc,tempMat,10);
    end;
%% % �������
Time_end = cputime;
Time_used = Time_end - Time_begin;    % ����������̺ķѵ�ʱ��
Time_clock = fix(clock);
disp(sprintf('�������ʱ�䣺 %d-%d-%d  %d:%d:%d',...
              Time_clock(1),Time_clock(2),Time_clock(3), ...
              Time_clock(4),Time_clock(5),Time_clock(6)));
disp(sprintf('���������������ʱ�䣺 Time_used = %g �롣',Time_used));
disp(sprintf('------�������------'));
disp(sprintf('  '));
disp(sprintf('  '));


%1.���ʵ����㷨�ܼ򵥣���֪����ʸ����ֱ����MVDR׼�����
%2.����ά���˲��Ǿ�������ʸ����������źţ�����������н��������������źŵĻ����ʸ��
%   ����ά�ɽ��������Ȩʸ����
%3.�����༶ά���˲��㷨�����ھ���ά���˲������Ͼ��н��ȴ���

