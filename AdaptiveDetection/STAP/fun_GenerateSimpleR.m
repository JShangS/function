function R=fun_GenerateSimpleR(K,N,Num,CNR)
% 2010.6.4
% ��ά��
% K��������
% N����Ԫ��
% Num: �Ƕȵķ���
% CNR:�����
% ����ʱ����Ӧ�źŴ������Ӳ�Э�������Ĳ��� 

% Num=181;
theta=linspace(-1/2,1/2,Num);
Beta=1;
Ac=(10^(CNR/10))^0.5; % ����������Ϊ1
Rc=zeros(K*N);
for i=1:length(theta)
    a=exp(1i*2*pi*(0:N-1)'*theta(i));
    b=exp(1i*2*pi*(0:K-1)'*Beta*theta(i));
    v=kron(b,a);
    Rc=Rc+v*v';
end
Rc=K*N*Rc/sum(eig(Rc))*Ac^2;
R=Rc+eye(K*N); 