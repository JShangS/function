function [ R,opta] = fun_KLPW(R0,R1,Rcut)%,
%FUN_KLPW 此处显示有关此函数的摘要
%   此处显示详细说明
%%R0:先验协防差
%%R1:SCM
%%R1:CUT单元协方差
%%%基于KL度量的最白化\alpha因子计算，参考《思考，4月8日》
[N,~] = size(R0);
% Rcut = fun_Positive(Rcut);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opta = min(1,abs(trace((Rcut - R1) / (R0 - R1)) / N));
% R = opta * R0 + (1 - opta) * R1;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 10;
R = Rcut;
for i = 1:iter
    opta = min(1,abs(trace((R - R1) / (R0 - R1)) / N));
    R = opta * R0 + (1 - opta) * R1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=0;
% iter_num = 10;
% for i = 1:iter_num
%     t1 = det((a*R0 + (1-a)*R1) \ Rcut);
%     a =  abs(trace((Rcut*t1 - R1)/(R0 - R1))/N);
% end
% opta = min(1,a);
% R = opta * R0 + (1 - opta) * R1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% da = 0.01;
% opta = 0;
% I = eye(N);
% mind = 1000000;
% for a = da:da:1
%     R_t =  a * R0 + (1-a) * R1;
% %     t1 = R_t^(-0.5) * Rcut * R_t^(-0.5);
%     t1 = R_t^(-1) * Rcut ;
%     min_det_t = trace((t1) - I) -log(det(t1));
%     if mind > min_det_t
%         mind = min_det_t;
%         opta = a;
%         R = R_t;
%     end
% end
end

