function [ R_rho ] = fun_rho( rho,N)
%FUN_ROU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%�������Э��������
%%%rho����������
%%%N������ʸ��ά��
R_rho = zeros(N,N);
for i=1:N
    for j=1:N
        R_rho(i,j)=rho^abs(i-j);
    end
end
end

