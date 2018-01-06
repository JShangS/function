function [ alpha,beta ] = fun_IG_ML(r)
%%《Maximum Likelihood Estimation forCompound-Gaussian Clutter with Inverse GammaTexture》
%%(15) (14) 
%%复合高斯invert Gamma 纹理参数ML估计
% alpha:形状参数
% beta：尺度参数
Ns = length(r);
beta_t = 1.:0.001:2; %%%IPIX
r2 = abs(r).^2;
L_beta_min = 1e10;
for i = 1:length(beta_t)
    t1 = sum(r2./(beta_t(i)*r2+1));
    t2 = sum(log(beta_t(i)*r2+1));
    L_beta = (Ns*beta_t(i)*t1)/(Ns-beta_t(i)*t1)-t2;
    L_beta = abs(L_beta);
    if L_beta<L_beta_min
        L_beta_min = L_beta;
        beta = beta_t(i);
    end
end
t3 = sum(r2./(beta*r2+1));
alpha = Ns/(beta*t3)-1;
end

