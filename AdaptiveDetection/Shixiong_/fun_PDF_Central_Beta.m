function PDF=fun_PDF_Central_Beta(N,M,Beta)
% 刘维建.2012.10.19
% 用到自编函数 fun_IG.m
% N、M : the number of degrees of freedom (DOF's)
% delt2: the noncentral parameter
% Beta can be a vector or a scalar.
% delt2 can be a vector, but it must have the same length with F.

% This script is written according to equation (29) in "Performance of the
% adaptive sidelobe blanker detection algorithm in homogeneous environments" C.D. Richmond, 2000
% or equation (A2-12) in "Adaptive detection and parameter estimation for
% multidimensional signal models", Kelly,1989

% PDF=factorial(N+M-1)/(factorial(N-1)*factorial(M-1))*Beta.^(N-1).*(1-Beta).^(M-1);
% PDF=exp(log(factorial(N+M-1))-log(factorial(N-1))-log(factorial(M-1)))*Beta.^(N-1).*(1-Beta).^(M-1);
% 当N或M太大时，溢出，例如阶乘最大计算到 (170)!

PDF=1/beta(N,M)*Beta.^(N-1).*(1-Beta).^(M-1);

