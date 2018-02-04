clc
clear 
close all
open('IPIX.fig');
h_line=get(gca,'Children');%get line handles
xdata=get(h_line,'Xdata');
ydata=get(h_line,'Ydata');