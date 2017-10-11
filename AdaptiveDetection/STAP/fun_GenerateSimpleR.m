function R=fun_GenerateSimpleR(K,N,Num,CNR)
% 2010.6.4
% 刘维建
% K：脉冲数
% N：阵元数
% Num: 角度的份数
% CNR:杂噪比
% 《空时自适应信号处理》，杂波协方差矩阵的产生 

% Num=181;
theta=linspace(-1/2,1/2,Num);
Beta=1;
Ac=(10^(CNR/10))^0.5; % 设噪声功率为1
Rc=zeros(K*N);
for i=1:length(theta)
    a=exp(1i*2*pi*(0:N-1)'*theta(i));
    b=exp(1i*2*pi*(0:K-1)'*Beta*theta(i));
    v=kron(b,a);
    Rc=Rc+v*v';
end
Rc=K*N*Rc/sum(eig(Rc))*Ac^2;
R=Rc+eye(K*N); 