function[R_MKA]=fun_MKA(X,R_KA,x0)
%FUN_MKA此处显示有关此函数的摘要
%此处显示详细说明
%%%利用AIWCM/AML对先验协方差进行修正,AML无法修正
%%%
[M,L] = size(X);
%%AML修正
% tao_child = 1;
% tao_parent = 0;
% R_MKA = R_KA;
% while (abs(tao_child-tao_parent)>0.1)%
%     tao_parent = tao_child;
%     R_MKA_inv = inv(R_MKA);
%     tao_child = diag(abs(X'*R_MKA_inv*X)/M);
%     R_MKA_t = 0;
%     for i = 1:L
%         R_MKA_t = R_MKA_t+X(:,i)*X(:,i)'/L/tao_child(i);
%     end
%     R_MKA = abs(R_MKA_t);
% end
%%AIWCM修正
beta_child=ones(L,1);
beta_parent=zeros(L,1);
R_MKA=R_KA;
count = 1;
while(norm(beta_parent-beta_child,'fro')/norm(beta_child,'fro')>0.1)
    count;
    beta_parent=beta_child;
    R_KA_M_inv=inv(R_MKA);
    beta_child_t1=X'*R_KA_M_inv*x0;
    beta_child_t2=diag(X'*R_KA_M_inv*X);
    beta_child=abs(beta_child_t1./beta_child_t2);
    R_KA_M_t=0;
    for i=1:L
        R_KA_M_t=R_KA_M_t+beta_child(i)*(X(:,i)*X(:,i)')/sum(beta_child);%sum(beta_child);
    end
    R_MKA=(R_KA_M_t);
    count=count+1;
    if count>50
       break;
    end
end
%%%%
end

