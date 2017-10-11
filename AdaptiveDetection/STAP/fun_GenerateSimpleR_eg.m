clear all; close all; clc


Na=4;
Np=6;
N=Na*Np;
CNR=60;
Num=361;
R=fun_GenerateSimpleR(Na,Np,Num,CNR);
Eig_R=eig(R);
figure;plot(db(Eig_R),'b-x')
iR=inv(R);
Num=128;
beta=1;
theta=(-1/2+1/Num:1/Num:1/2).';
fd=beta*theta;

for i=1:length(theta)
    a=exp(j*2*pi*(0:Na-1).'*theta(i));
    for jj=1:length(fd)        
        b=exp(j*2*pi*(0:Np-1).'*fd(jj));
        steer=kron(b,a);   % steer=steer/norm(steer);
        res_opt(jj,i)=1/(steer'*iR*steer); 
    end
end
[Theta Fd]=meshgrid(theta,fd);
res_opt=res_opt/max(max(abs(res_opt)));
res_opt=10*log10(abs(res_opt));
% res_opt(find(res_opt<-120))=-120;
figure;mesh(Theta,Fd,res_opt)
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')
figure;contour(Theta,Fd,res_opt)
xlabel('方位');ylabel('Doppler');zlabel('功率 (dB)')

