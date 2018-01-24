function [ Train ] = fun_TrainData_IGCC( N,L,M,lambda,mu,opt_train)
%FUN_TRAINDATA_IGCC 此处显示有关此函数的摘要
%   此处显示详细说明
%JerryShang，2017.11.16
%%逆gamma纹理的复合高斯杂波
%%N,L行列,导向矢量维度，参考单元数目
%%M，复高斯协方差
%%lamda：形状参数
%%mu：尺度参数
%%opt:1: 每个距离单元的纹理值不一样，在CPI内每个距离单元的纹理值一样，SIRP
%%%%%%2：每个单元纹理值一样，为部分均匀环境
%%%%%%3：每个单元纹理值每时每刻都不一样，最广义的compound Gaussian
%%%详见刘军大论文
if nargin == 5
    opt_train = 1;
end
if opt_train == 1
    pct=sqrt(ones(L,1)./gamrnd(lambda,mu,L,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,1)*pct.').*X;
elseif opt_train == 2
    pct=sqrt(1./gamrnd(lambda,mu,1,1));
    X = fun_TrainData_gauss(N,L,M);
    Train=(ones(N,L)*pct.').*X;
elseif opt_train == 3
    pct = zeros(N,L);
    for i = 1:L
        pct(:,i)=sqrt(1./gamrnd(lambda,mu,N,1));
    end  
    X = fun_TrainData_gauss(N,L,M);
    Train=pct.*X;
else
    error('只有1，2, 3三种选择')
end

end

