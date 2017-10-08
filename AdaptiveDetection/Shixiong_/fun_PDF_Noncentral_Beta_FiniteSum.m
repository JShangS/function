function PDF=fun_PDF_Noncentral_Beta_FiniteSum(N,M,delta2,Beta)
% 刘维建.2016.06.07
% N、M : the number of degrees of freedom (DOF's)
% delta2: the noncentral parameter
% Beta can be a vector or a scalar.
% delta2 can be a vector, but it must have the same length with F.

% This script is written according  "Adaptive detection and parameter estimation for
% multidimensional signal models", Kelly,1989

% PDF=factorial(N+M-1)/(factorial(N-1)*factorial(M-1))*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delt2).*hypergeom(N+M,M,delt2.*(1-Beta));
% PDF=1/beta(N,M)*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delt2).*hypergeom(N+M,M,delt2.*(1-Beta));
% PDF=1/beta(N,M)*Beta.^(N-1).*(1-Beta).^(M-1).*exp(-delta2).*hypergeom(N+M,M,delta2.*(1-Beta));
% 当N或M太大时，溢出，例如阶乘最大计算到 (170)!
tmp=zeros(size(Beta));
for k=0:N
    num=delta2.^k.*(1-Beta).^(M+k-1);
    den=factorial(k)*factorial(N-k)*factorial(M+k-1);
    tmp=tmp+num/den;
end
% PDF=N*factorial(N+M-1)*exp(-delta2.*Beta).*Beta.^(N-1).*tmp;
lnPDF=log(N)+log(factorial(N+M-1))+(-delta2.*Beta)+log(Beta.^(N-1))+log(tmp);
PDF=exp(lnPDF); %【注】：当非中心参数太大，例如40dB时，就会出现e^(NNN)计算，而Matlab能够计算的最大量纲小于e^710，e^709是可以计算出来的
% figure; semilogy( lnPDF )
% 上式是把Kelly1989年报告中的(A2-12)式代入(A2-23)式得到的

