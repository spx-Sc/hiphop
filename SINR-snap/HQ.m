function [ r_temp ] = HQ( R,M,L )
%MDL 此处显示有关此函数的摘要
%   此处显示详细说明
[V,Deig]=eig(R);
%         for i=1:M
%             value_D(1,i) = Deig(i,i); % 将特征值排成一行
%         end;
        value_D=diag(Deig).';
        SortValueD = sort(value_D,'descend');
        for k=1:M-1
            DiagSum(k) = sum(SortValueD(k+1:M))/(M-k);
            DiagProd(k) = (prod(SortValueD(k+1:M)))^(1/(M-k));
            alpha(k) = DiagSum(k)/DiagProd(k);
            HQ(k) = (M-k)*L*log(alpha(k))+k/2*(2*M-k)*log(log(L));
        end;
        [d_min,r_temp] = min(HQ);

end

