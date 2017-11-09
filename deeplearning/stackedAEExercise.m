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

%��һ���Զ��������w1��b1
stack{1}.w = reshape(sae1OptTheta(1:hiddenSizeL1*inputSize), ...
                     hiddenSizeL1, inputSize);
stack{1}.b = sae1OptTheta(2*hiddenSizeL1*inputSize+1:2*hiddenSizeL1*inputSize+hiddenSizeL1);
%�ڶ����Զ��������w2��b2
stack{2}.w = reshape(sae2OptTheta(1:hiddenSizeL2*hiddenSizeL1), ...
                     hiddenSizeL2, hiddenSizeL1);
stack{2}.b = sae2OptTheta(2*hiddenSizeL2*hiddenSizeL1+1:2*hiddenSizeL2*hiddenSizeL1+hiddenSizeL2);

% Initialize the parameters for the deep model
% [params, netconfig] = stack2params(stack)
% �����ǽ�stack��ε���������������Ƕ��������ת����һ������params��
% ����������ʹ�ø����Ż��㷨�������Ż�������
% Netconfig�б�����Ǹ�����������Ϣ��
% ����netconfig.inputsize��ʾ��������������ڵ�ĸ�����
% netconfig.layersizes�е�Ԫ�طֱ��ʾÿһ���������Ӧ�ڵ�ĸ�����

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
% �ú����ڲ�ʵ������������ʧ��������ʧ������ÿ������ƫ���ļ��㡣
% ������ʧ�����Ǹ�ʵ��ֵ����Ȼ��ֻ��1���ˣ�����㷽���Ǹ���sofmax������������ģ�ֻ��֪����ǩֵ��softmax������ֵ���ɡ�
% ����ʧ���������в�����ƫ��ȴ�кܶ�������ÿ��������Ӧ�þ���һ��ƫ��ֵ����Щ�������������˶��������ģ����һ�������softmax�Ǹ������ġ�
% ����softmax�ǲ��ֵ�ƫ���Ǹ����乫ʽֱ�ӻ�ã������������ǲ�����ͨ��BP�㷨��������õ������ȼ���ÿһ������ֵ��Ȼ�����ø����ֵ�������w��b����
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

%!!!!!!!!!!!!!�������ǳ���Ҫ������ʶ����ֻ��%91��
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
% �����������ʵ���Ƕ������data���ݽ���Ԥ�⣬����data��Ӧ���������Ƕ��١�
%   ����thetaΪ��������Ĳ����������˷��������ֵ����磩��numClassesΪ�����������netconfigΪ����Ľṹ������

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
