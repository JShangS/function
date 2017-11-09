clc
clear
close all
load trainDataWarhead;
load trainDataOther;

trainData=zeros(784,60000);
for a=1:392
    trainData(a,:)=trainDataWarhead(a,:);
    trainData(a+392,:)=trainDataOther(a,:);
end
trainLabels=zeros(60000,1);
for i=1:30000
    trainLabels(i,1)=2;
    trainLabels(i+30000,1)=1;
end
trainData2d = zeros(28,28,60000);
for i = 1:28
    trainData2d(1:end,i,1:end)=trainData((i-1)*28+1:i*28,1:end);
end
trainData2d = mat2gray(trainData2d,[0,255]);
trainData2d = fftshift(trainData2d,1);
imagesc(trainData2d(:,:,1))
save('JS_trainData2d.mat','trainData2d','trainLabels')
% mesh(fftshift(trainData2d(:,:,1),1))
% view(0,-90)