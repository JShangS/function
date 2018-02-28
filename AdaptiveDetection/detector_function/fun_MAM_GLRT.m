function [ T_MAM_GLRT ] = fun_MAM_GLRT( MAM,x0,p,lambda,mu )
%  ���� <Knowledge-based adaptive detection of radar targets in generalized Pareto clutter>
%%%�ļ����������MAMģ�ͣ�Ȼ��GLRT���������Gamma�����µ������Ȼ����
%%%%MAM:��ģ��
%%%p:Ŀ�굼��ʸ��
%%%x0:CUT
%%%lambda����gamma�����������״������
%%%mu����gamma����������߶Ȳ�����
if nargin<4 %%%%Ĭ�ϲ���Ϊ3,1
    lambda = 3;
    mu = 1;
end
[N,~,L]=size(MAM); 
t1 = zeros(L,1);%%%1������ϵ��
t0 = zeros(L,1);%%%0������ϵ��
for i=1:L
    MAM(:,:,i) = inv(MAM(:,:,i));
    a=(p'*MAM(:,:,i)*x0)/(p'*MAM(:,:,i)*p);
    t1(i) = N/(lambda*mu*(x0-a*p)'*MAM(:,:,i)*(x0-a*p));
    t0(i) = N/(lambda*mu*x0'*MAM(:,:,i)*x0);
end
%%Э�������
%%1�����µ�
t1 = reshape(t1,1,1,L);
t1 = repmat(t1,N,N,1);
R_MAM1 = sum(t1.*MAM,3)/L;
%%0�����µ�
t0 = reshape(t0,1,1,L);
t0 = repmat(t0,N,N,1);
R_MAM0 = sum(t0.*MAM,3)/L;
%%%�����
a = (p'*R_MAM1*x0)/(p'*R_MAM1*p);
f1 = det(R_MAM1)*((x0-a*p)'*R_MAM1*(x0-a*p)+1/mu)^(-lambda-N);
f0 = det(R_MAM0)*(x0'*R_MAM0*x0+1/mu)^(-lambda-N);
T_MAM_GLRT = abs(f1/f0);
end

