%% CS294A/CS294W Stacked Autoencoder Exercise

%  Instructions
%  ------------
% 
%  This file contains code that helps you get started on the
%  sstacked autoencoder exercise. You will need to complete code in
%  stackedAECost.m
%  You will also need to have implemented sparseAutoencoderCost.m and 
%  softmaxCost.m from previous exercises. You will need the initializeParameters.m
%  loadMNISTImages.m, and loadMNISTLabels.m files from previous exercises.
%  
%  For the purpose of completing the assignment, you do not need to
%  change the code in this file. 
%
%%======================================================================
%% STEP 0: Here we provide the relevant parameters values that will
%  allow your sparse autoencoder to get good filters; you do not need to 
%  change the parameters below.
clc; clear ; close all;


inputSize = 28 * 28;
% numClasses = 10;
numClasses = 2;
hiddenSizeL1 = 200;    % Layer 1 Hidden Size
hiddenSizeL2 = 200;    % Layer 2 Hidden Size
sparsityParam = 0.1;   % desired average activation of the hidden units.
                       % (This was denoted by the Greek alphabet rho, which looks like a lower-case "p",
		               %  in the lecture notes). 
lambda = 3e-3;         % weight decay parameter       
beta = 3;              % weight of sparsity penalty term       

%%======================================================================
%% STEP 1: Load data from the MNIST database
%
%  This loads our training data from the MNIST database files.

% Load MNIST database files
% trainData = loadMNISTImages('train-images.idx3-ubyte');
% trainData2 = loadMNISTImages('tiain-images.idx3-ubyte');
% trainData=zeros(784,60000);
% for a=1:3000
%     trainData(:,a)=trainData1(:,a);
%     trainData(:,a+30000)=trainData2(:,a+30000);
% end
% trainLabels = loadMNISTLabels('train-labels.idx1-ubyte');
% load trainDataWarhead;
% load trainDataOther;
% 
% trainData=zeros(784,60000);
% for a=1:392
%     trainData(a,:)=trainDataWarhead(a,:);
%     trainData(a+392,:)=trainDataOther(a,:);
% end
% trainLabels=zeros(60000,1);
% for i=1:30000
%     trainLabels(i,1)=2;
%     trainLabels(i+30000,1)=1;
% end
load
% trainLabels(trainLabels == 0) = 10; % Remap 0 to 10 since our labels need to start from 1

%%======================================================================
%% STEP 2: Train the first sparse autoencoder
%  This trains the first sparse autoencoder on the unlabelled STL training
%  images.
%  If you've correctly implemented sparseAutoencoderCost.m, you don't need
%  to change anything here.


%  Randomly initialize the parameters
sae1Theta = initializeParameters(hiddenSizeL1, inputSize);

%% ---------------------- YOUR CODE HERE  ---------------------------------
%  Instructions: Train the first layer sparse autoencoder, this layer has
%                an hidden size of "hiddenSizeL1"
%                You should store the optimal parameters in sae1OptTheta


%  Use minFunc to minimize the function
addpath minFunc/
options.Method = 'lbfgs'; % Here, we use L-BFGS to optimize our cost
                          % function. Generally, for minFunc to work, you
                          % need a function pointer with two outputs: the
                          % function value and the gradient. In our problem,
                          % sparseAutoencoderCost.m satisfies this.
options.maxIter = 200;	  % Maximum number of iterations of L-BFGS to run 
options.display = 'on';
options.Corr = 10;

[sae1OptTheta, cost] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                   inputSize, hiddenSizeL1, ...
                                   lambda, sparsityParam, ...
                                   beta, trainData), ...
                                   sae1Theta, options);


save 'sae1OptTheta.mat' sae1OptTheta

load sae1OptTheta;








% -------------------------------------------------------------------------



%%======================================================================
%% STEP 2: Train the second sparse autoencoder
%  This trains the second sparse autoencoder on the first autoencoder
%  featurse.
%  If you've correctly implemented sparseAutoencoderCost.m, you don't need
%  to change anything here.

[sae1Features] = feedForwardAutoencoder(sae1OptTheta, hiddenSizeL1, ...
                                        inputSize, trainData);

%  Randomly initialize the parameters
sae2Theta = initializeParameters(hiddenSizeL2, hiddenSizeL1);

%% ---------------------- YOUR CODE HERE  ---------------------------------
%  Instructions: Train the second layer sparse autoencoder, this layer has
%                an hidden size of "hiddenSizeL2" and an inputsize of
%                "hiddenSizeL1"
%
%                You should store the optimal parameters in sae2OptTheta




%  Use minFunc to minimize the function
addpath minFunc/
options.Method = 'lbfgs'; % Here, we use L-BFGS to optimize our cost
                          % function. Generally, for minFunc to work, you
                          % need a function pointer with two outputs: the
                          % function value and the gradient. In our problem,
                          % sparseAutoencoderCost.m satisfies this.
options.maxIter = 200;	  % Maximum number of iterations of L-BFGS to run 
options.display = 'on';


[sae2OptTheta, cost] = minFunc( @(p) sparseAutoencoderCost(p, ...
                                   hiddenSizeL1, hiddenSizeL2, ...
                                   lambda, sparsityParam, ...
                                   beta, sae1Features), ...
                                   sae2Theta, options);
save 'sae2OptTheta.mat' sae2OptTheta
load sae2OptTheta;







% -------------------------------------------------------------------------


%%======================================================================
%% STEP 3: Train the softmax classifier
%  This trains the sparse autoencoder on the second autoencoder features.
%  If you've correctly implemented softmaxCost.m, you don't need
%  to change anything here.

[sae2Features] = feedForwardAutoencoder(sae2OptTheta, hiddenSizeL2, ...
                                        hiddenSizeL1, sae1Features);

%  Randomly initialize the parameters
saeSoftmaxTheta = 0.005 * randn(hiddenSizeL2 * numClasses, 1);


%% ---------------------- YOUR CODE HERE  ---------------------------------
%  Instructions: Train the softmax classifier, the classifier takes in
%                input of dimension "hiddenSizeL2" corresponding to the
%                hidden layer size of the 2nd layer.
%
%                You should store the optimal parameters in saeSoftmaxOptTheta 
%
%  NOTE: If you used softmaxTrain to complete this part of the exercise,
%        set saeSoftmaxOptTheta = softmaxModel.optTheta(:);


options.maxIter = 100;
softmaxModel = softmaxTrain(hiddenSizeL2, numClasses, lambda, ...
                            sae2Features, trainLabels, options);


saeSoftmaxOptTheta = softmaxModel.optTheta(:)
save 'saeSoftmaxOptTheta.mat' saeSoftmaxOptTheta
load saeSoftmaxOptTheta;

% -------------------------------------------------------------------------



%%======================================================================
%% STEP 5: Finetune softmax model

% Implement the stackedAECost to give the combined cost of the whole model
% then run this cell.

% Initialize the stack using the parameters learned
stack = cell(2,1);

%第一层自动编码机的w1，b1
stack{1}.w = reshape(sae1OptTheta(1:hiddenSizeL1*inputSize), ...
                     hiddenSizeL1, inputSize);
stack{1}.b = sae1OptTheta(2*hiddenSizeL1*inputSize+1:2*hiddenSizeL1*inputSize+hiddenSizeL1);
%第二层自动编码机的w2，b2
stack{2}.w = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), ...
                     hiddenSizeL2, hiddenSizeL1);
stack{2}.b = sae2OptTheta(2*hiddenSizeL2*hiddenSizeL1+1:2*hiddenSizeL2*hiddenSizeL1+hiddenSizeL2);

% Initialize the parameters for the deep model
% [params, netconfig] = stack2params(stack)
% 　　是将stack层次的网络参数（可能是多个参数）转换成一个向量params，
% 这样有利用使用各种优化算法来进行优化操作。
% Netconfig中保存的是该网络的相关信息，
% 其中netconfig.inputsize表示的是网络的输入层节点的个数。
% netconfig.layersizes中的元素分别表示每一个隐含层对应节点的个数。

[stackparams, netconfig] = stack2params(stack);
stackedAETheta = [ saeSoftmaxOptTheta ; stackparams ];

%% ---------------------- YOUR CODE HERE  ---------------------------------
%  Instructions: Train the deep network, hidden size here refers to the '
%                dimension of the input to the classifier, which corresponds 
%                to "hiddenSizeL2".
%
%

DEBUG = false;

if DEBUG
    lambda = 0;
    hiddenSizeL2 = 2;
    trainData = trainData(:,1:10);
    trainLabels = trainLabels(1:10);
    netconfig.layersizes = {};
    netconfig.layersizes = [netconfig.layersizes;2];
    netconfig.layersizes = [netconfig.layersizes;2];
    stackedAETheta = stackedAETheta(1:1596);
end

% [cost, grad] = stackedAECost(stackedAETheta, inputSize, hiddenSizeL2, ...
%                                               numClasses, netconfig, ...
%                                               lambda, trainData, trainLabels);


if DEBUG
%     numGrad = computeNumericalGradient( @(x) stackedAECost(x, numClasses, ...
%                                     inputSize, lambda, trainData, trainLabels), stackedAETheta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ cost, grad ] = stackedAECost(theta, inputSize, hiddenSize, numClasses, netconfig,lambda, data, labels)
% 该函数内部实现整个网络损失函数和损失函数对每个参数偏导的计算。
% 其中损失函数是个实数值，当然就只有1个了，其计算方法是根据sofmax分类器来计算的，只需知道标签值和softmax输出层的值即可。
% 而损失函数对所有参数的偏导却有很多个，因此每个参数处应该就有一个偏导值，这些参数不仅包括了多个隐含层的，而且还包括了softmax那个网络层的。
% 其中softmax那部分的偏导是根据其公式直接获得，而深度网络层那部分这通过BP算法方向推理得到（即先计算每一层的误差值，然后利用该误差值计算参数w和b）。
    numGrad = computeNumericalGradient( @(x) stackedAECost(x, inputSize, hiddenSizeL2, ...
                                              numClasses, netconfig, lambda, ...
                                              trainData, trainLabels), stackedAETheta);
                                
    % Use this to visually compare the gradients side by side
    disp([numGrad grad]); 

    % Compare numerically computed gradients with those computed analytically
    diff = norm(numGrad-grad)/norm(numGrad+grad);
    disp(diff); 
    % The difference should be small. 
    % In our implementation, these values are usually less than 1e-7.

    % When your gradients are correct, congratulations!
end

%!!!!!!!!!!!!!这个步骤非常重要，否则识别率只有%91！
lambda = 1e-4;

addpath minFunc/
options.Method = 'lbfgs'; 
options.maxIter = 400;	 
options.display = 'on';

[stackedAEOptTheta, cost] = minFunc( @(p) stackedAECost(p, inputSize, hiddenSizeL2, ...
                                              numClasses, netconfig, lambda, ...
                                              trainData, trainLabels), stackedAETheta, options);
save 'stackedAEOptTheta.mat' stackedAEOptTheta;

load stackedAEOptTheta;

% -------------------------------------------------------------------------



%%======================================================================
%% STEP 6: Test 
%  Instructions: You will need to complete the code in stackedAEPredict.m
%                before running this part of the code
%

% Get labelled test images
% Note that we apply the same kind of preprocessing as the training set
testData = loadMNISTImages('t10k-images.idx3-ubyte');
% testLabels = loadMNISTLabels('t10k-labels.idx1-ubyte');
% 
% testLabels(testLabels == 0) = 10; % Remap 0 to 10
testLabels=zeros(10000,1);
for i=1:10000
    testLabels(i,1)=2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pred] = stackedAEPredict(theta, inputSize, hiddenSize, numClasses, netconfig, data)
% 　这个函数其实就是对输入的data数据进行预测，看该data对应的输出类别是多少。
%   其中theta为整个网络的参数（包括了分类器部分的网络），numClasses为所需分类的类别，netconfig为网络的结构参数。

[pred] = stackedAEPredict(stackedAETheta, inputSize, hiddenSizeL2, ...
                          numClasses, netconfig, testData);

acc = mean(testLabels(:) == pred(:));
fprintf('Before Finetuning Test Accuracy: %0.3f%%\n', acc * 100);

[pred] = stackedAEPredict(stackedAEOptTheta, inputSize, hiddenSizeL2, ...
                          numClasses, netconfig, testData);

acc = mean(testLabels(:) == pred(:));
fprintf('After Finetuning Test Accuracy: %0.3f%%\n', acc * 100);

% Accuracy is the proportion of correctly classified images
% The results for our implementation were:
%
% Before Finetuning Test Accuracy: 87.7%
% After Finetuning Test Accuracy:  97.6%
%
% If your values are too low (accuracy less than 95%), you should check 
% your code for errors, and make sure you are training on the 
% entire data set of 60000 28x28 training images 
% (unless you modified the loading code, this should be the case)
