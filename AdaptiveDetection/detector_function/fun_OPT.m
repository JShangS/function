function [ Topt ] = fun_OPT( R,x0,p )
%%%���ż������Э���������֪�����
%%%R��Covariance Matrix��x0��CUT��p��steering vector
iR = inv(R);
Topt = abs(p'*iR*x0)^2 / (p'*iR*p) / (x0'*iR*x0);
end

