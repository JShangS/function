% Time Doppler plot
function [logTD,time,doppl]=tdoppl(nc,I,Q,rangebin,xpix,ypix)
global file_chosen;
if nargin<5, xpix=250; end
if nargin<6, ypix=80; end
if rangebin>0
    NRangeBins = 1;%只循环一次
else%对所有距离单元遍历
    NRangeBins = size(I,1);
end

figure('Name','TD','Position',[0,0,1280,960]);
title('Time Doppler');

range=netcdf.getVar(nc,netcdf.inqVarID(nc,'range'));

for range_cell_index=1:NRangeBins;
    image_index=mod(range_cell_index,14);
    if(image_index==0) image_index=14; end
    Ir=I(range_cell_index,:);
    Qr=Q(range_cell_index,:);
    % load I and Q data
    %[I,Q]=ipixcdf(nc,txpol,rangebin);
    N=size(I,2);
    
    % adjust fft window and averaging to match image size
    wdw=max(128,2^ceil(log(4*N/xpix)/log(2)));
    yavg=max(1,ceil(wdw/ypix));
    timeStep=wdw/4;
    M=floor((N-wdw)/timeStep);
    
    % short time fourier transforms
    hw=hamming(wdw);
    TD=zeros(floor(wdw/yavg),M);
    y=zeros(yavg,floor(wdw/yavg)); leny=prod(size(y));
    for m=1:M,
        i=(m-1)*timeStep;
        x=Ir(1+i:wdw+i)+j*Qr(1+i:wdw+i);
        x=hw'.*x;
        fftx=abs(fftshift(fft(x)));
        y(:)=fftx(1:leny);
        TD(:,m)=mean(y,1)';
    end
    
    % time and normalized frequency
    PRF = netcdf.getVar(nc,netcdf.inqVarID(nc,'PRF'));;%nc{'PRF'}(1);              % Pulse Repetition Frequency [Hz]
    % time=[0 M]/PRF;
    % freq=[-0.5 0.5]*PRF;
    time=(wdw/2+(0:M-1)*timeStep)/PRF;
    freq=((0:wdw/yavg-1)/(wdw/yavg)-0.5)*PRF;
    
    % convert to doppler velocity
    RF_Freq= netcdf.getVar(nc,netcdf.inqVarID(nc,'RF_frequency'));
    doppl=freq*3e8/(2*RF_Freq*1e9);
    
    % enhance log-plot with noise floor
    mn=mean(TD,2);mn_list(range_cell_index)=mean(mn);
    [mn,indx]=sort(mn);
    noise=TD(indx(1:2),:);
    noiseFloor=median(noise(:));
    TD(find(TD<noiseFloor))=noiseFloor;
    logTD=log(TD);
    mn=min(min(logTD)); mx=max(max(logTD));
    logTD=(mn_list(range_cell_index)/mn_list(1))*(logTD-mn)*63/(mx-mn);
    % display image
    % if nargout<1,
%     subplot(5,3,image_index);
    image(time,doppl,logTD); set(gca,'ydir','normal');
    title(['Range ' num2str(range_cell_index) ': ' num2str( range(range_cell_index)) ' m']);
    % end
    
%     if image_index==14 || range_cell_index==NRangeBins
% %         subplot(5,3,15);
%         plot(mn_list,'o:');
%         title('Mean of Log Intensity');xlabel('Rangebin');ylabel('Mean');
% %         pause;
%     end
end

