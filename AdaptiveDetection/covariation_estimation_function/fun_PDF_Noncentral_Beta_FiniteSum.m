function PDF=fun_PDF_Noncentral_Beta_FiniteSum(N,M,delta2,Beta)
% ��ά��.2016.06.07
% N��M : the number of degrees of freedom (DOF's)
% delta2: the noncentral parameter
% Beta can be a vector or a scalar.
% delta2 can be a vector, but it must have the same length with F.

% This script is written according  "Adaptive detection and parameter estimation for
% multidimensional signal models", Kelly,1989

% PDF=factorial(N+M-1)/(factorial(N-1)*factorial(M-1))*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delt2).*hypergeom(N+M,M,delt2.*(1-Beta));
% PDF=1/beta(N,M)*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delt2).*hypergeom(N+M,M,delt2.*(1-Beta));
% PDF=1/beta(N,M)*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delta2).*hypergeom(N+M,M,delta2.*(1-Beta));
% ��N��M̫��ʱ�����������׳������㵽 (170)!
tmp=zeros(size(Beta));
for k=0:N
    num=delta2.^k.*(1-Beta).^(M+k-1);
    den=factorial(k)*factorial(N-k)*factorial(M+k-1);
    tmp=tmp+num/den;
end
% PDF=N*factorial(N+M-1)*exp(-delta2.*Beta).*Beta.^(N-1).*tmp;
lnPDF=log(N)+log(factorial(N+M-1))+(-delta2.*Beta)+log(Beta.^(N-1))+log(tmp);
PDF=exp(lnPDF); %��ע�����������Ĳ���̫������40dBʱ���ͻ����e^(NNN)���㣬��Matlab�ܹ�������������С��e^710��e^709�ǿ��Լ��������
% figure; semilogy( lnPDF )
% ��ʽ�ǰ�Kelly1989�걨���е�(A2-12)ʽ����(A2-23)ʽ�õ���

