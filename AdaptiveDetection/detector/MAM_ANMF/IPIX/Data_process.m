Global_varity;
load(matFile) 
N = 8;
Range = 14;
[row,col] = size(sig);
if col>row
    sig = sig.';
end
[M,L] = size(sig);
% figure()
% imagesc(abs(sig))
sig_test_t = sig(2:N+1,:);
[X,Y]=meshgrid(1:L,linspace(-0.5,0.5,M));
MTD = abs(fftshift(fft(sig,[],1)));
% figure()
% mesh(X,Y,MTD);
R_KA = zeros(N,N);
L_KA = 7500;
for i = 1:L_KA
    X_KA = sig((i-1)*N+1:i*N,Range);
    R_KA = R_KA+X_KA*X_KA'/L_KA;
end
% figure(3)
% mesh(abs(R_KA))
save(matFile,'R_KA','Range','N','sig');
