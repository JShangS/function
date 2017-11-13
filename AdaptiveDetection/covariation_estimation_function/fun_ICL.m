function [ R_ICL ] = fun_ICL( p,z,RKA,R,mu,lamda,H)
%FUN_ICL 此处显示有关此函数的摘要
%   此处显示详细说明
%%RKA：先验协方差
%%R:训练数据协方差
%%mu，lamda，逆gamma分布的参数
%%H：1，H1下的估计结果，H0下的估计结果
%%p:导向矢量，z:CUT，待检测单元
[N,~] = size(z);
delt_min= 100;
beta_real = 0;
for beta = 0.01:0.01:1
    a = (p'*inv(beta*(RKA-R)+R)*z)/(p'*inv(beta*(RKA-R)+R)*p);
    if H==1
       t1 = log(det(beta*(RKA-R)+R));
       t2 = -(lamda+N)*log((z-a*p)'*(beta*(RKA-R)+R)*(z-a*p)+1/mu);
       delt = t1 + t2;
    else
       t1 = log(det(beta*(RKA-R)+R));
       t2 = -(lamda+N)*log(z'*(beta*(RKA-R)+R)*z+1/mu);
       delt = t1 + t2;
    end
    if abs(delt)<delt_min
        delt_min = abs(delt);
        beta_real = beta;
    end
end
 R_ICL = abs(beta_real*inv(RKA) + (1-beta_real)*inv(R));


% 
% if H == 1
%     
%     a_1 = p'*RKA*z/(p'*RKA*p);
%     a_2 = p'*R*z/(p'*R*p); 
%     b1 = (z-a_1*p)'*RKA*(z-a_1*p);
%     b2 = (z-a_2*p)'*R*(z-a_2*p);
%     alpha = -N/(mu*b1*(N+(lamda+N)*det(RKA)));
%     beta = -N/(mu*b2*(N+(lamda+N)*det(R)));
% else
%     b1 = z'*RKA*z;
%     b2 = z'*R*z;
%     alpha = -N/(mu*b1*(N+(lamda+N)*det(RKA)));
%     beta = -N/(mu*b2*(N+(lamda+N)*det(R)));
% end
% alpha = abs(alpha);
% beta = abs(beta);
% sum_all = abs(alpha+beta);
% alpha = abs(alpha/sum_all)
% beta = abs(beta/sum_all)
% R_ICL = abs(alpha*inv(RKA) + beta*inv(R));




% beta_child = 0;
% beta_parent = 1; 
% a = (p'*R_ICL_inv*z)/(p'*R_ICL_inv*p);
% count = 1;
% 
% if H==1
%    fai = abs(((z-a*p)'*R_inv*(z-a*p)+1/mu)/((z-a*p)'*(RKA_inv-R_inv)*(z-a*p)));
% else
%    fai = abs((z'*R_inv*z+1/mu)/(z'*(RKA_inv-R_inv)*z));
% end
% delt_min = 1000;
% beta_real = 0;
% for beta = 0.01:0.01:1
%     t1 = log(abs(det(RKA_inv-R_inv)/det(beta*(RKA_inv-R_inv)+R_inv)))
%     t2 = -(lamda+N)*log((beta+fai)^(-1));
%     delt = t1+t2;
%     if abs(delt)<delt_min
%         delt_min = delt;
%         beta_real = beta
%     end
% end
% R_ICL_inv = beta_real*RKA_inv + (1-beta_real)*R_inv;


% while abs(beta_child-beta_parent)>0.1
%     beta_parent = beta_child
%     a = (p'*R_ICL*z)/(p'*R_ICL*p);
%     if H==1
%         fai = abs(((z-a*p)'*R*(z-a*p)+1/mu)/((z-a*p)'*(RKA-R)*(z-a*p)));
%     else
%         fai = abs((z'*R*z+1/mu)/(z'*(RKA-R)*z));
%     end
%     beta_child = min(abs(((lamda+N)*fai*det(RKA-R))/det(beta_parent*(RKA-R)+R)),1)
%     R_ICL = beta_child*RKA + (1-beta_child)*R;
%     count = count +1;
%     if count >10
%         break;
%     end
% end
end

