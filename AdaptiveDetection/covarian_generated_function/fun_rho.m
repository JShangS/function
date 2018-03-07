function [ R_rho ] = fun_rho( rho,N,fd)
%FUN_ROU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%�������Э��������
%%%rho����������
%%%N������ʸ��ά��
R_rho = zeros(N,N);
if nargin<3
    fd = 0;
end
L = length(rho);
for l = 1:L
    for i=1:N
        for j=1:N
            R_rho(i,j,l)=rho(l)^abs(i-j)*exp(1j*2*pi*fd*(i-j));
        end
    end
end
end


