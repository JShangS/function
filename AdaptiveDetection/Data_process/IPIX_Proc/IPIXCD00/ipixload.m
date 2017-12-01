% function [I,Q,meanIQ,stdIQ,inbal]=ipixload(nc,pol,rangebin,mode);
%
% Loads I and Q data from IPIX radar cdf file.
% Inputs:
%   nc       - pointer to netCDF file
%   pol      - Transmit/receive polarization ('vv','hh','vh', or 'hv')
%   rangebin - Range bin(s) to load
%   mode     - Pre-processing mode:
%              'raw' does not apply any corrections to the data;
%              'auto' applies automatic corrections, assumes that radar
%                 does not look at still objects, such as land;
%              'dartmouth' first removes land, using knowledge of the geometry
%                 of the Dartmouth site. Then same as 'auto'.
%
% Outputs:
%   I,Q      - Pre-processed in-phase and quadrature component of data
%   meanIQ   - Mean of I and Q used in pre-processing
%   stdIQ    - Standard deviation of I and Q used in pre-processing
%   inbal    - Phase inbalance [degrees] used in pre-processing
%
% Rembrandt Bakker, McMaster University, 2001-2002
%

function [I,Q,meanIQ,stdIQ,inbal,adc_data]=ipixload(nc,pol,rangebin,mode)

%% check inputs
varRange = netcdf.getVar(nc,netcdf.inqVarID(nc,'range'));
nrange=length(varRange);
if rangebin==0,
    rangebin=1:nrange;
end

if rangebin<1 | rangebin>nrange,
    disp(['Warning: rangebin ' num2str(rangebin) ...
        ' not found in file ' nc.NetCDF_file_name(1:end-1)]);
    return;
end

%% in some cdf files, the unsigned flag is not set correctly %%
adc_data= netcdf.getVar(nc,netcdf.inqVarID(nc,'adc_data'));

H_txpol     = 1;
V_txpol     = 2;
% Like_adc_I  = nc{'adc_like_I'}(1)+1;
% Like_adc_Q  = nc{'adc_like_Q'}(1)+1;
% Cross_adc_I = nc{'adc_cross_I'}(1)+1;
% Cross_adc_Q = nc{'adc_cross_Q'}(1)+1;
Like_adc_I  = 1;%netcdf.getVar(nc,netcdf.inqVarID(nc,'adc_like_I'))+1;
Like_adc_Q  = 2;%netcdf.getVar(nc,netcdf.inqVarID(nc,'adc_like_Q'))+1;
Cross_adc_I = 3;%netcdf.getVar(nc,netcdf.inqVarID(nc,'adc_cross_I'))+1;
Cross_adc_Q = 4;%netcdf.getVar(nc,netcdf.inqVarID(nc,'adc_cross_Q'))+1;
%% extract correct polarization from cdffile %%
pol=lower(pol);
if length(size(adc_data))==3,
    % read global attribute TX_polarization
    txpol=netcdf.getAtt(nc,-1,'TX_polarization');  %nc.TX_polarization(1);
    if pol(1)~=lower(txpol),
        fname=strrep(netcdf.getAtt(nc,-1,'NetCDF_file_name'),'\0','');%nc.NetCDF_file_name(:)
        disp(['Warning: file ' fname ' does not contain ''' ...
            txpol(1) ''' transmit polarization.']);
    end
    switch pol,
        case {'hh','vv'}, xiq=adc_data([Like_adc_I Like_adc_Q],rangebin,:);
        case {'hv','vh'}, xiq=adc_data([Cross_adc_I Cross_adc_Q],rangebin,:);
    end
    I=double(permute(xiq(1,:,:),[ 2 3 1]));
    Q=double(permute(xiq(2,:,:),[2 3 1 ]));
    %   figure(1)
    %   subplot(2,2,1);
    %   plot(I(1,:));
    %   subplot(2,2,2);
    %   plot(Q(1,:));
    %
else
    switch pol,
        case 'hh', xiq=adc_data([1 3],rangebin,H_txpol,:);
        case 'hv', xiq=adc_data([1 2],rangebin,H_txpol,:);
        case 'vv', xiq=adc_data([2 4],rangebin,V_txpol,:);
        case 'vh', xiq=adc_data([2 3],rangebin,V_txpol,:);
    end
    I=double(permute(xiq(1,:,1,:),[2 4 1 3]));
    Q=double(permute(xiq(2,:,1,:),[2 4 1 3]));
end

%% apply corrections to I and Q data %%
if nargin<3, mode='auto'; end
switch mode,
    case 'raw',
        meanIQ=[0 0]; stdIQ=[1 1]; inbal=0;
    case 'auto',
        % Pre-processing
        meanIQ=mean([I(:) Q(:)]); 
        stdIQ=std([I(:) Q(:)],1);
        I=(I-meanIQ(1))/stdIQ(1);
        Q=(Q-meanIQ(2))/stdIQ(2);
        sin_inbal=mean(I(:).*Q(:));
        inbal=asin(sin_inbal)*180/pi;
        I=(I-Q*sin_inbal)/sqrt(1-sin_inbal^2);
    case 'dartmouth',
        % Define rectangular patches of land in Dartmouth campaign.
        % Format: [azmStart azmEnd  rangeStart rangeEnd]
        landcoord=[
            0  70     0  600;
            305 360     0  600;
            30  55     0 8000;
            210 305     0 4700;
            320 325  2200 2700;
            ];
        % Exclude land from data used to estimate pre-processing parameters
        azm=mod(ipixazm(nc),360);
        range=varRange(rangebin);
        nbin=length(rangebin);
        ok=ones(size(I));
        for i=1:size(landcoord,1),
            for r=1:nbin,
                if range(r)>=landcoord(i,3) & range(r)<=landcoord(i,4),
                    ok(find(azm>=landcoord(i,1) & azm<=landcoord(i,2)),r)=0;
                end
            end
        end
        ok=find(ok);
        if length(ok)<100,
            disp(['Warning: not enough sweeps for land-free pre-processing.']);
            ok=ones(size(I));
        end
        % Pre-processing
        meanIQ=mean([I(ok) Q(ok)]); stdIQ=std([I(ok) Q(ok)],1);
        I=(I-meanIQ(1))/stdIQ(1);
        Q=(Q-meanIQ(2))/stdIQ(2);
        sin_inbal=mean(I(ok).*Q(ok));
        inbal=asin(sin_inbal)*180/pi;
        I=(I-Q*sin_inbal)/sqrt(1-sin_inbal^2);
        
end

%imagesc(abs(double(I(:,:))+j*double(Q(:,:))));
[nPol nrangebins nsweep]=size(adc_data);
range=varRange(:);

% %%  PPI.m
% rpix=200;
% mode=''; 
% seaonly=strcmp(mode,'seaonly');
% 
% % averaging to match window size
% 
% xavg=max(1,ceil(nsweep/(4*pi*rpix)));
% nrange=length(range);
% maxrange=max(range);
% 
% % process radar returns
% AA=zeros(nrange,floor(nsweep/xavg));
% ampl=zeros(xavg,floor(nsweep/xavg)); lenAmpl=prod(size(ampl));
% for r=1:nrange,
%     absiq=abs(I+j*Q)*sqrt(prod(stdIQ));
%     ampl(:)=absiq(1:lenAmpl);
%     AA(r,:)=mean(ampl.^2,1);
%     if seaonly,
%         % correct for power drop off
%         AA(r,:)=AA(r,:)*(range(r)/range(end))^3;
%     end
% end
% 
% % convert to polar coordinates
% az=zeros(xavg,floor(nsweep/xavg)); lenAz=prod(size(az));
% azan=netcdf.getVar(nc,netcdf.inqVarID(nc,'azimuth_angle'));
% az(:)=unwrap(azan(1:lenAz)*pi/180);
% az=mean(az,1);
% % remove corrupt azimuth angles
% daz=median(diff(az));
% indx=find(sign(diff(az))~=sign(daz));
% if ~isempty(indx),
%     az=az(1:indx(1)); AA=AA(:,1:indx(1));
% end
% if seaonly,
%     % remove land
%     indx=find(az>70*pi/180 & az<215*pi/180);
%     az=az(indx); AA=AA(:,indx);
% end
% ppr=rpix/maxrange;
% xPolar=ppr*range*cos(0.5*pi-az);
% yPolar=ppr*range*sin(0.5*pi-az);
% [xI,yI]=meshgrid(min(xPolar(:))-8:max(xPolar(:))+8,...
%     min(yPolar(:))-8:max(yPolar(:))+8);
% xI=round(xI); yI=round(yI);
% rI=sqrt(xI.^2+yI.^2);
% azI=0.5*pi-angle(xI+j*yI); indx=find(azI<min(az)); azI(indx)=azI(indx)+2*pi;
% ok=find(rI>=ppr*min(range) & rI<=ppr*maxrange & azI>=min(az) & azI<=max(az));
% AAI=interp2(az,ppr*range,AA,azI(ok),rI(ok));
% logA=zeros(size(xI)); logA(ok)=log(AAI);
% 
% % add rticks
% indx=find(rI>=2 & rI<=4);
% logA(indx)=1;
% rtick=[0:10:350]*pi/180;
% indx=find(rtick>min(az) & rtick<max(az));
% rtick=rtick(indx);
% for r=1:length(rtick)
%     indx1=find((rI>rpix & rI<=rpix+8));
%     indx2=find(azI(indx1)>rtick(r)-0.75./rI(indx1) & azI(indx1)<rtick(r)+0.75./rI(indx1));
%     logA(indx1(indx2))=1;
% end
% rtick=[0:90:350]*pi/180;
% indx=find(rtick>min(az) & rtick<max(az));
% rtick=rtick(indx);
% for r=1:length(rtick)
%     indx1=find(((rI>rpix & rI<=rpix+8) | (rI>4 & rI<ppr*min(range))));
%     indx2=find(azI(indx1)>rtick(r)-1.5./rI(indx1) & azI(indx1)<rtick(r)+1.5./rI(indx1));
%     logA(indx1(indx2))=1;
% end
% 
% if seaonly,
%     % scale for best visual contrast
%     ea=360-mean(nc{'elevation_angle'}(:));
%     rmin=30/((ea+0.45)*pi/180);
%     rmax=30/(ea*pi/180);
%     best=find(rI>=rpix*rmin/maxrange & rI<=rpix*rmax/maxrange & ...
%         azI>=min(az) & azI<=max(az));
%     if length(best)<0.01*length(ok), best=ok; end
%     sortA=sort(logA(best));
%     len=length(sortA);
%     med=sortA(round(len/2));
%     mad=median(abs(sortA-med));
%     mn=med-3*mad; mx=med+3*mad;
% else
%     mn=min(logA(ok)); mx=max(logA(ok));
% end
% indx=find(logA(ok)<mn); logA(ok(indx))=mn;
% indx=find(logA(ok)>mx); logA(ok(indx))=mx;
% logA(ok)=2+(logA(ok)-mn)*61/(mx-mn);
% 
% if seaonly, map=gray; map(:,3)=map(:,3)*0.9; else map=jet; end
% map(1,:)=[1 1 1];
% map(2,:)=[0 0 0];
% 
% 
%     colormap(map);
%     image(logA+1,'cdatamapping','direct'); set(gca,'ydir','normal');
%     set(gcf,'units','pixels','position',[60 20 100+size(logA,2) 80+size(logA,1)]);
%     set(gca,'units','pixels','position',[60 60 size(logA,2) size(logA,1)]);
% 
%     