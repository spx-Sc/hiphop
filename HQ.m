function [ r_temp ] = HQ( R,M,L )
%MDL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[V,Deig]=eig(R);
%         for i=1:M
%             value_D(1,i) = Deig(i,i); % ������ֵ�ų�һ��
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

