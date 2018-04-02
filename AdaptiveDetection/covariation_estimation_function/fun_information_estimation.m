function [ R_MAM, distance,ratio] = fun_information_estimation(R0, MAM,opt,alpha)
%FUN_INFORMATION_ESTIMATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%%%%%%������Ϣ���ζ�����Э������Ʒ�����
%%%R0:�ο�Э����
%%%MAM:��ģ�����Э����ΪN*N��LΪģ�͸���
%%%����ѡ�� r��ReimanDistance��l: Log Euclidean Distance
%%%e: Euclidean Distance, ro:Root Euclidean Distance
%%%c: CholeskyDistance, p: Power-Euclidean distance

if nargin==2
    opt = 'r';%ReimanDistance
    alpha = 1;
elseif nargin == 3
    alpha = 1;
end
[N,~,L]=size(MAM); 
distance = zeros(L,1);
switch opt
    case 'r'%ReimanDistance
        for i=1:L
            distance(i) = fun_ReimanDistance(R0,MAM(:,:,i));
        end
%         distance = distance/sum(distance);
    case 'l' %Log Euclidean Distance
        for i=1:L
            distance(i) = fun_LogED(R0,MAM(:,:,i));
        end
%         distance = distance/sum(distance);
    case 'c' %CholeskyDistance
        for i=1:L
            distance(i) = fun_CholeskyDistance(R0,MAM(:,:,i));
        end
%         distance = distance/sum(distance);
    case 'e'%Euclidean Distance
        for i=1:L
            distance(i) = fun_EuclideanDistance(R0,MAM(:,:,i));
        end
%         distance = distance/sum(distance);
    case 'p'%Power-Euclidean distance
        for i=1:L
            distance(i) = fun_PowerED(R0,MAM(:,:,i),alpha);
        end
    case 'ro'%Root Euclidean Distance
        for i=1:L
            distance(i) = fun_RootED(R0,MAM(:,:,i));
        end
    case 'kl' %KL divergence
        for i=1:L
            distance(i) = fun_KLD(R0,MAM(:,:,i));
        end
    case 'bh' %Bhattacharyya distance
        for i=1:L
            distance(i) = fun_BhD(R0,MAM(:,:,i));
        end   
   case 'h' %Hellinger distance
        for i=1:L
            distance(i) = fun_HD(R0,MAM(:,:,i));
        end 
   case 'skl' %sKL divergence
        for i=1:L
            distance(i) = fun_sKLD(R0,MAM(:,:,i));
        end
        
    otherwise
end
h = 1;
% distance = distance/sum(distance);
Sum_ratio = sum(exp(-(distance)/h^2));
ratio = exp(-(distance))/Sum_ratio;
ratio = reshape(ratio,1,1,L);
ratio = repmat(ratio,N,N,1);
R_MAM = sum(ratio.*MAM,3);
ratio = ratio(1,1,1:L);
end

