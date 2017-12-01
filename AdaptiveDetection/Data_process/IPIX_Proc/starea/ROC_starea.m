clear all
close all
clc

% data1 主目标21
% ncid=netcdf.open('19980303_195104_ANTSTEP.cdf');
% data2 主目标22
ncid=netcdf.open('19931107_135603_starea.cdf');

nrange=14;
nsweep=131072;
sp=50;
pre=1e-5;
step=1;
step2=1e1;
N=nsweep/sp;
y=zeros(nrange,nsweep);
yabs=zeros(nrange,nsweep);

for bin=1:nrange
    [ I , Q ]=ipixload(ncid , 'hh' , bin , 'auto');
    y(bin,:)=conj(I+1i*Q)';
end

yest=y(9,1:sp:end);
yestabs=abs(yest);

ydec=y(10,1:sp:end);
ydecabs=abs(ydec);

[ k , beta ] = ComARCHEst( yestabs , pre , step );
[ p0 , p1 , sig2 , delta ] = ComNARCHEst( yestabs , pre , step2 );

yt12=[mean(yestabs)  yestabs(1:end-1)].^2;

v_gau=var(real(yest));
v_arch=0.5*(k+beta*yt12);
v_narch =0.5*( p0 * sig2^delta + p1 * yt12.^delta).^( 1 / delta) ;

p_f=logspace(-3,0);
n=length(p_f);
PD=zeros(3,n);
N2=length(ydecabs);

for i=1:n
    
    lamda_gau=sqrt(-v_gau*log(p_f(i)));
    lamda_arch=sqrt(-v_arch*log(p_f(i)));
    lamda_narch=sqrt(-v_narch*log(p_f(i)));
    
    n_gau=sum((ydecabs-lamda_gau)>0);
    n_arch=sum((ydecabs-lamda_arch)>0);
    n_narch=sum((ydecabs-lamda_narch)>0);
    
    PD(1,i)=n_gau/N2;
    PD(2,i)=n_arch/N2;
    PD(3,i)=n_narch/N2;
    
end
figure
semilogx(p_f,PD(1,:) , 'Color' , 'blue' , 'Marker' , 'o');
hold on
semilogx(p_f,PD(2,:) , 'Color' , 'red' , 'Marker' , '+');
hold on
semilogx(p_f,PD(3,:) , 'Color' , 'green' , 'Marker' , '*' );
legd=legend('Gaussian','ARCH','NARCH','Location','NorthWest');
xlabel('Probability of False Alarm');
ylabel('Probability of Detection')
ylim([0 1])
grid on