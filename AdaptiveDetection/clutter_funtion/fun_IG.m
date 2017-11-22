function IG=fun_IG(m,a)
% ¡ıŒ¨Ω®£¨2011.8.9
% incomplete Gamma function difined in
% "Performance of the Adaptive Sidelobe Blanker Detection Algorithm in Homogeneous Enviroments"
% m is a scalar, a can be a vector.

% temp=0;
% for k=0:m-1
%     temp=temp+a^k/factorial(k);
% end
% IG=exp(-a)*temp;
temp=zeros(size(a));
for k=0:m-1
    temp=temp+a.^k/factorial(k);
end
IG=exp(-a).*temp;
% figure;plot(fun_IG(5,1:100))
% figure;plot(fun_IG(5,0:0.01:1))