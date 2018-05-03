function [ R_NS ] = fun_NSCMC( X,opt )
%X:ѵ������
%%��һ������Э�������
%һ����һ�����뵥Ԫ
if nargin==1
    opt=1;
end
[M,N] = size(X);
NX = zeros(M,N);
R_NS = zeros(M,M);
if opt==1
    %%Adaptive matched filter detection in spherically invariant noise
    %%���ݹ�һ��
    for i = 1:N
        NX(:,i) = X(:,i)/sqrt(norm(X(:,i),'fro')^2/M);
        R_NS = R_NS + fun_Positive(NX(:,i),4)/N;
    end
    
elseif opt==2
    %%Performance analysis of two covariance matrix estimators in
    %%compound-Gaussian clutter
    for i = 1:N
         t = fun_Positive(X(:,i),4)/(X(:,i)'*X(:,i))/M;
        R_NS = R_NS + t/N;
    end   
end
end

