function [ Tamf ] = fun_AMF( R,x0,p )
%FUN_AMF �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%AMF
%%%R��Э���x0��CUT��p������ʸ��
iR = inv(R);
Tamf = abs(p'*iR*x0)^2/abs(p'*iR*p);   
end

