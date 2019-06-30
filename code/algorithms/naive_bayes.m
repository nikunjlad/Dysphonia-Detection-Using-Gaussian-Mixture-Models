%% naive bayes

load('data2.mat');
load('finalNB.mat');

% SVMModel = fitcnb(Xtr, ytr,'DistributionNames', 'kernel', ...
%     'OptimizeHyperparameters','all',...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus'));

% SVMModel = fitcnb(Xtr, ytr,'DistributionNames','kernel','Width',0.00059773,'Kernel','box');


%% Performance parameters for training validating data

% fprintf('\n--------------- Training Data ---------------\n');
% [labeltr,~] = predict(SVMModel,Xtr);
% results(labeltr,ytr,'Training');
%   
% %% Performance parameters for cross validating data
% 
% fprintf('\n--------------- Cross Validation Data ---------------\n');
% [labelcv,~] = predict(SVMModel,Xcv);
% results(labelcv,ycv,'CV');
%   
% %% Performance parameters for testing data
% 
% fprintf('\n--------------- Testing Data ---------------\n');
% labeltst = predict(SVMModel,Xtst);
% results(labeltst,ytst,'Testing');

%% Performance parameters for training entire dataset

fprintf('\n--------------- Entire Dataset ---------------\n');
[labeltt,~] = predict(SVMModel,X);
results(labeltt,y,'Training');