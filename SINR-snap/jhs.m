function aa=jhs(w,Rxp,u,bi,ybsn,C)
%using the method of J2!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    result=0;
    for i=1:120
        result=result+fi(i)*(u(i)^2);
    end
    aa=0.5*w'*Rxp*w+0.5*result+bi;



%     result=0;%using the method of J1!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         for i=1:120
%             if (abs(u(i))<ybsn)%using the method of J1!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 adder=0;
%             else
%                 adder=abs(u(i)-ybsn);
%             end
%             result=result+adder;
%         end
% 
%     aa=0.5*w'*Rxp*w+C*result;
end