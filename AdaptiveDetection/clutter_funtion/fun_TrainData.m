function [ X ] = fun_TrainData( varargin )
%FUN_TRAINDATA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%JerryShang��2017.11.17
%%%ʵ�ֺ������أ������ɲ�ͬ�ķֲ���ѵ������
%%'k':K�ֲ���
%%'g'��gauss,
%%'p'��generalized Pareto clutter,����gamma�����ĸ��ϸ�˹�ֲ�
str = varargin{1}; %%%�ж�Ҫʵ�ֵķֲ�����:
switch str
    case 'g'
        if nargin < 4
           error('Gauss�ֲ�ѵ�������������Ϊ4����g ѡ�����ʸ��ά����ѵ�����ݳ��ȣ�Э����'); 
        end
        N = varargin{2};
        L = varargin{3};
        R = varargin{4}; %%Э����
        X = fun_TrainData_gauss(N,L, R);
    case 'k'
        if nargin < 5
           error('K�ֲ�ѵ�������������Ϊ5����k ѡ�����ʸ��ά����ѵ�����ݳ��ȣ�Э����,��״����'); 
        end
        N = varargin{2};
        L = varargin{3};
        R = varargin{4};
        v = varargin{5};
        X = fun_TrainData_K( N,L, R, v);
    case 'p'
        if nargin < 6
           error('IG�ֲ�ѵ�����������������6����p ѡ�����ʸ��ά����ѵ�����ݳ��ȣ�Э����,��״����,�߶Ȳ�����SIRPѡ��'); 
        end
        N = varargin{2};
        L = varargin{3};
        R = varargin{4};
        lamda = varargin{5};
        mu = varargin{6};
        if nargin == 6
            opt_train = 1;
            X = fun_TrainData_IGCC( N,L,R,lamda,mu,opt_train);
        else
           opt_train = varargin{7};
           X = fun_TrainData_IGCC( N,L,R,lamda,mu,opt_train);
        end
        
    otherwise
            error('����ֻ��3��');
end
end
