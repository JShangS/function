function [ R_CC,alpha0 ] = fun_CCIter(X,R_KA )
%FUN_CC �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%������⣬ѵ���������Ƶ�Э���������Э�����������ϣ�����͹�Ż��õ����ϵ����
%%X:ѵ������
%R,�������Ƶ�Э����
%R_KA:����Э����
R = fun_SCMN(X);
[R_CC,alpha0_1] = fun_CC(X,R_KA);
while(1)
    [R_CC,alpha0] = fun_CC2(X,R_CC,R_KA);
    if abs(alpha0_1 - alpha0) <1e-5
        break;
    end
    alpha0_1 = alpha0;
end
% iter = 20;
% alpha0 = zeros(iter+1,1);
% [R_CC,alpha0(1)] = fun_CC(X,R,R_KA);
% for i = 1 :iter
%     [R_CC,alpha0(i+1)] = fun_CC(X,R_CC,R_KA);
% end
end


