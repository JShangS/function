clear
close
clc


ncid=netcdf.open('19931107_135603_starea.cdf');
tdoppl(ncid,'hh',9,250,80);


[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

%各变量的值
for i=0:nvars-1
    assignin('base',netcdf.inqVar(ncid,i),netcdf.getVar(ncid,i));
end

xpix=3e2;
ypix=1e2;
PRF=1e3;

% averaging to match window size
nsweep=size(adc_data,4);
xavg=max(1,ceil(nsweep/xpix));
nrange=length(range);

% process radar returns
A=zeros(nrange,floor(nsweep/xavg));
ampl=zeros(xavg,floor(nsweep/xavg));
lenAmpl=numel(ampl);
time=(1:size(A,2))/PRF*xavg;

polar={'hh','hv','vv','vh'};
for i =1:length(polar)
    for r=1:nrange,
        [ I , Q , mnIQ , sdIQ , inbal ] = ipixload( ncid , char(polar(i)) , r , 'auto' );
        absiq = abs( I + 1i*Q ) * sqrt( prod( sdIQ ) );
        ampl(:)=absiq(1:lenAmpl);
        A(r,:)=mean(ampl,1);
    end
    
    % range averaging
    while size(A,1)>ypix,
        A=0.5*(A(1:2:end-1,:)+A(2:2:end,:));
        range=0.5*(range(1:2:end-1)+range(2:2:end));
    end
    
    % indexed color log plot
    logA=log(A);
    mn=min(min(logA));
    mx=max(max(logA));
    logA=(logA-mn)*63/(mx-mn);
    subplot(4,1,i);
    image(time,range,logA);
    set(gca,'ydir','normal');
    title(char(polar(i)));
end

