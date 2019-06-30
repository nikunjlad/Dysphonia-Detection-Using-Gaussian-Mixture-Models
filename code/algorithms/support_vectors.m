%% SVM without i-vectors

%% Loading the dataset and the trained RBF Model
load('dt12.mat');
load('md1.mat');

% [Xtr,ytr,Xcv,ycv,Xtst,ytst,X,y] = randomize(X,y);
%% SVM Classifier Model (NOTE: Uncomment it only if you are training it for the first time, after which saved model loaded)

% SVMModel = fitcsvm(Xtr, ytr,'Standardize',true,'KernelFunction','rbf',...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus'));

% SVMModel = fitcsvm(Xtr, ytr,'Standardize',true,'KernelFunction','rbf',...
%     'BoxConstraint',990,'KernelScale',1.2378);
%  
%% Performance parameters for training validating data

% fprintf('\n--------------- Training Data ---------------\n');
% [labeltr,~] = predict(SVMModel,Xtr);
% [~,~,~,~,~,~,~,~,~,~] = results(labeltr,ytr,'Training',0);
  
% %% Performance parameters for cross validating data
% 
% fprintf('\n--------------- Cross Validation Data ---------------\n');
% [labelcv,~] = predict(SVMModel,Xcv);
% [~,~,~,~,~,~,~,~,~,~] = results(labelcv,ycv,'CV',0);
%   
% %% Performance parameters for testing data
% 
% fprintf('\n--------------- Testing Data ---------------\n');
% labeltst = predict(SVMModel,Xtst);
% [~,~,~,~,~,~,~,~,~,~] = results(labeltst,ytst,'Testing',0);

%% Performance parameters for training entire dataset

fprintf('\n--------------- Entire Dataset ---------------\n');
[labeltt,~] = predict(SVMModel,X);
[~,~,~,~,~,~,~,~,~,~] = results(labeltt,y,'Entire Dataset',1);






