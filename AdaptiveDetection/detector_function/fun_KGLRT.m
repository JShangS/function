function [ Tglrt ] = fun_KGLRT( R,x0,p )
%FUN_KGLRT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%Kelly��1S-GLRT
%%%R��Э���x0��CUT��p������ʸ��
Tamf =  fun_AMF(R,x0,p);
iR = inv(R);
tmp=abs(x0'*iR*x0);
Tglrt = Tamf/(1+tmp);
end

