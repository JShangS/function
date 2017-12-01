% function [azm,outlier]=ipixazm(nc)
%
% In the Dartmouth database there are a few 
% files in which the azimuth angle
% time-series is corrupt, in particular
% data sets 8, 16, 23, 24, 28, 110, 254, 276.
% This function fixes this bug.
% Inputs:
%   nc       - pointer to netCDF file
%
% Outputs:
%   azm      - corrected azimuth angle time series. 
%   outliers - indices of azimuth angles considered outliers
%
% Rembrandt Bakker, McMaster University, 2001-2002
%

function [azm,outlier]=ipixazm(nc)

azm=unwrap(nc{'azimuth_angle'}(:)*pi/180);
dazm=diff(azm);
meddazm=median(dazm);
outlier=find(abs(dazm)>2*abs(meddazm));
dazm(outlier)=meddazm;
azm=cumsum([azm(1); dazm])*180/pi;


