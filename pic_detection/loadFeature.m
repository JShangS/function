clc
clear 
close all
%%直接loade训练好的feature
load cnnPooledFeatures.mat
load MTDData_fu20dB.mat
trainImages = trainData;
[~,~,numImages_JS] = size(trainImages);
numTrainImages = numImages_JS;
numTestImages =  numImages_JS;
%% STEP 4: Use pooled features for classification
%  Now, you will use your pooled features to train a softmax classifier,
%  using softmaxTrain from the softmax exercise.
%  Training the softmax classifer for 1000 iterations should take less than
%  10 minutes.

% Add the path to your softmax solution, if necessary
% addpath /path/to/solution/

% Setup parameters for softmax
softmaxLambda = 1e-4;
% numClasses = 4;
numClasses = 2;%--JS--
% trainLabels = trainLabels;%--JS--
% Reshape the pooledFeatures to form an input vector for softmax
softmaxX = permute(pooledFeaturesTrain, [1 3 4 2]);
softmaxX = reshape(softmaxX, numel(pooledFeaturesTrain) / numTrainImages,...
    numTrainImages);
softmaxY = trainLabels;

options = struct;
options.maxIter = 200;
softmaxModel = softmaxTrain(numel(pooledFeaturesTrain) / numTrainImages,...
    numClasses, softmaxLambda, softmaxX, softmaxY, options);

%%======================================================================
%% STEP 5: Test classifer
%  Now you will test your trained classifer against the test images
% testLabels = testLabels;%--JS--
softmaxX = permute(pooledFeaturesTest, [1 3 4 2]);
softmaxX = reshape(softmaxX, numel(pooledFeaturesTest) / numTestImages, numTestImages);
softmaxY = testLabels;

[pred] = softmaxPredict(softmaxModel, softmaxX);
acc = (pred(:) == softmaxY(:));
acc = sum(acc) / size(acc, 1);
fprintf('Accuracy: %2.3f%%\n', acc * 100);

% You should expect to get an accuracy of around 80% on the test images.