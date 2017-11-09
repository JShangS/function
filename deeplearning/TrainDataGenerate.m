% clc;
% clear all;
% testData = loadMNISTImages('t10k-images.idx3-ubyte');
% testLabels = loadMNISTLabels('t10k-labels.idx1-ubyte');
% trainData = loadMNISTImages('train-images.idx3-ubyte');
% trainLabels = loadMNISTLabels('train-labels.idx1-ubyte');
% trainDataWarhead=zeros(392,60000);
% for i=1:392
%     trainDataWarhead(i,:)=trainData(i+181,:);
% end
% save 'trainDataOther.mat' trainDataOther;
clc
clear all
close all

load qiusxh_5F;
load qiusxh_10F;
load qiusxh_20F;
load qiusxhF;
load qiusxv_5F;
load qiusxv_10F;
load qiusxv_20F;
load qiusxvF;
load sc1h_5F;
load sc1hF;
load sc1v_5F;
load sc1vF;
load sxh_5F;
load sxhF;
load sxv_5F;
load sxvF;
load zhuisxh_5F;
load zhuisxh_10F;
load zhuisxh_20F;
load zhuisxhF;
load zhuisxv_5F;
load zhuisxv_10F;
load zhuisxv_20F;
load zhuisxvF;
trainData=zeros(900,4424);
trainData(1:900,1:101)=sc1h_5F(1:900,:);
trainData(1:900,102:202)=sc1hF(1:900,:);
trainData(1:900,203:303)=sc1v_5F(1:900,:);
trainData(1:900,304:404)=sc1vF(1:900,:);
trainData(1:900,405:605)=sxh_5F(1:900,:);
trainData(1:900,606:806)=sxhF(1:900,:);
trainData(1:900,807:1007)=sxv_5F(1:900,:);
trainData(1:900,1008:1208)=sxvF(1:900,:);
trainData(1:900,1209:1409)=qiusxh_5F(1:900,:);
trainData(1:900,1410:1610)=qiusxh_10F(1:900,:);
trainData(1:900,1611:1811)=qiusxh_20F(1:900,:);
trainData(1:900,1812:2012)=qiusxhF(1:900,:);
trainData(1:900,2013:2213)=qiusxv_5F(1:900,:);
trainData(1:900,2214:2414)=qiusxv_10F(1:900,:);
trainData(1:900,2415:2615)=qiusxv_20F(1:900,:);
trainData(1:900,2616:2816)=qiusxvF(1:900,:);
trainData(1:900,2817:3017)=zhuisxh_5F(1:900,:);
trainData(1:900,3018:3218)=zhuisxh_10F(1:900,:);
trainData(1:900,3219:3419)=zhuisxh_20F(1:900,:);
trainData(1:900,3420:3620)=zhuisxhF(1:900,:);
trainData(1:900,3621:3821)=zhuisxv_5F(1:900,:);
trainData(1:900,3822:4022)=zhuisxv_10F(1:900,:);
trainData(1:900,4023:4223)=zhuisxv_20F(1:900,:);
trainData(1:900,4224:4424)=zhuisxvF(1:900,:);

load sc1h_10F;
load sc1v_10F;
load sxh_10F;
load sxv_10F;
testData=zeros(900,604);
testData(1:900,1:101)=sc1h_10F(1:900,:);
testData(1:900,102:202)=sc1v_10F(1:900,:);
testData(1:900,203:403)=sxh_10F(1:900,:);
testData(1:900,404:604)=sxv_10F(1:900,:);





    