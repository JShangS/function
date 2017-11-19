function [ R ] = fun_TrainData( varargin )
%FUN_TRAINDATA 此处显示有关此函数的摘要
%   此处显示详细说明
%JerryShang，2017.11.17
%%%实现函数重载，可生成不同的分布的训练样本
%%'k':K分布，
%%'g'：gauss,
%%'p'：generalized Pareto clutter,即逆gamma纹理的复合高斯分布
str = varargin{1}; %%%判断要实现的分布类型:
switch str
    case 'g'
        if nargin ~= 4
           error('Gauss分布训练数据输入参数为4个：g 选项，导向矢量维数，训练数据长度，协方差'); 
        end
        N = varargin{2};
        L = varargin{3};
        M = varargin{4};
        R = fun_TrainData_gauss(N,L, M);
    case 'k'
        if nargin ~= 5
           error('K分布训练数据输入参数为5个：k 选项，导向矢量维数，训练数据长度，协方差,形状参数'); 
        end
        N = varargin{2};
        L = varargin{3};
        M = varargin{4};
        v = varargin{5};
        R = fun_TrainData_K( N,L, M, v);
    case 'p'
        if nargin < 6
           error('IG分布训练数据输入参数至少6个：p 选项，导向矢量维数，训练数据长度，协方差,形状参数,尺度参数，SIRP选项'); 
        end
        N = varargin{2};
        L = varargin{3};
        M = varargin{4};
        lamda = varargin{5};
        mu = varargin{6};
        if nargin == 6
            opt = 1;
        end
        R = fun_TrainData_IGCC( N,L,M,lamda,mu,opt);
    otherwise
            error('现在只有3种');
end
end

