function [ cost, grad ] = stackedAECost(theta, inputSize, hiddenSize, ...
                                              numClasses, netconfig, ...
                                              lambda, data, labels)
                                         
% stackedAECost: Takes a trained softmaxTheta and a training data set with labels,
% and returns cost and gradient using a stacked autoencoder model. Used for
% finetuning.
                                         
% theta: trained weights from the autoencoder
% visibleSize: the number of input units
% hiddenSize:  the number of hidden units *at the 2nd layer*
% numClasses:  the number of categories
% netconfig:   the network configuration of the stack
% lambda:      the weight regularization penalty
% data: Our matrix containing the training data as columns.  So, data(:,i) is the i-th training example. 
% labels: A vector containing labels, where labels(i) is the label for the
% i-th training example


%% Unroll softmaxTheta parameter

% We first extract the part which compute the softmax gradient
softmaxTheta = reshape(theta(1:hiddenSize*numClasses), numClasses, hiddenSize);

% Extract out the "stack"
stack = params2stack(theta(hiddenSize*numClasses+1:end), netconfig);

% You will need to compute the following gradients
softmaxThetaGrad = zeros(size(softmaxTheta));
stackgrad = cell(size(stack));
for d = 1:numel(stack)
    stackgrad{d}.w = zeros(size(stack{d}.w));
    stackgrad{d}.b = zeros(size(stack{d}.b));
end

cost = 0; % You need to compute this

% You might find these variables useful
M = size(data, 2);
groundTruth = full(sparse(labels, 1:M, 1));


%% --------------------------- YOUR CODE HERE -----------------------------
%%几个分类就几个
z2 = bsxfun(@plus, stack{1}.w*data, stack{1}.b);
a2 = sigmoid(z2);
% z3 = bsxfun(@plus, stack{2}.w*a2, stack{2}.b);
% a3 = sigmoid(z3);
zz = bsxfun(@plus, stack{2}.w*a2, stack{2}.b);
aa = sigmoid(zz);
% zz = bsxfun(@plus, stack{3}.w*a3, stack{3}.b); 
% aa = sigmoid(zz);

z4 = softmaxTheta * aa;
a4 = exp(z4);
%可以用这种方法替换我前面写的repmat
a4 = bsxfun(@rdivide, a4, sum(a4));

delta4 = -(groundTruth - a4);

% deltazz=(softmaxTheta' * delta4) .* sigmoidGrad(zz);
delta3=(softmaxTheta' * delta4) .* sigmoidGrad(zz);

% delta3 = (stack{3}.w' * deltazz) .* sigmoidGrad(z3);


delta2 = (stack{2}.w' * delta3) .* sigmoidGrad(z2);

softmaxThetaGrad = -(1. / M) * (groundTruth - a4) * aa'  + lambda * softmaxTheta;

% stackgrad{3}.w = (1. / M) * deltazz * a3' + lambda * stack{3}.w;
% stackgrad{3}.b = (1. / M) * sum(deltazz, 2);
stackgrad{2}.w = (1. / M) * delta3 * a2' + lambda * stack{2}.w;
stackgrad{2}.b = (1. / M) * sum(delta3, 2);
stackgrad{1}.w = (1. / M) * delta2 * data' + lambda * stack{1}.w;
stackgrad{1}.b = (1. / M) * sum(delta2, 2);

% cost = (1. / M) * sum((1. / 2) * sum((a4 - groundTruth).^2));

cost = -(1. / M) * sum(sum(groundTruth .* log(a4))) + (lambda / 2.) * ...
    sum(sum(softmaxTheta.^2)) + (lambda / 2.) * (sum(sum(stack{1}.w .^2)) ...
    + sum(sum(stack{2}.w .^2)));%+sum(sum(stack{3}.w .^2))

function grad = softmaxGrad(x)
    e_x = exp(-x);
    grad = e_x ./ (1 + (1-e_x).*e_x ).^2;
end

function grad = sigmoidGrad(x)
    e_x = exp(-x);
    grad = e_x ./ ((1 + e_x).^2); 
end


% -------------------------------------------------------------------------

%% Roll gradient vector
grad = [softmaxThetaGrad(:) ; stack2params(stackgrad)];

end


% You might find this useful
function sigm = sigmoid(x)
    sigm = 1 ./ (1 + exp(-x));
end
