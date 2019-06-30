load('V.mat');
% load('t17.mat');
X = v;
l = lasso(v,y);
% loss = zeros(size(l,2)-2,2);

% for i = 1:size(l,2)-2
    idx = find(l(:,36));
    A = X(:,idx);

    [h,~] = size(X);
    p1 = ceil(0.7 * h);
    p2 = p1 + ceil(0.2 * h);

    Xtr = X(1:p2,:);
    ytr = y(1:p2,:);

    Xtst = X(p2:end,:);
    ytst = y(p2:end,:);

%% SVM Classifier Model (NOTE: Uncomment it only if you are training it for the first time, after which saved model loaded)

    SVMModel = fitcknn(Xtr, ytr, 'Standardize', true, 'KernelScale', 'auto', 'KernelFunction', 'rbf',...
        'CacheSize', 8000);

    CVSVMModel = crossval(SVMModel);

    loss = kfoldLoss(CVSVMModel,'mode','individual');
% end


H = CVSVMModel.Trained{2};
traindex = training(CVSVMModel.Partition,2);
XTrain = X(traindex,:);
YTrain = y(traindex,:);
[labeltr,~] = predict(SVMModel,XTrain);
acc = mean(double(labeltr == YTrain)) * 100;
fprintf('\nTraining Set Accuracy: %f\n', acc);

cvdex = test(CVSVMModel.Partition,2);
Xcv = X(cvdex,:);
Ycv = y(cvdex,:);
[labelcv,~] = predict(SVMModel,Xcv);
acc = mean(double(labelcv == Ycv)) * 100;
fprintf('\nTraining Set Accuracy: %f\n', acc);

% [labeltst,~] = predict(SVMModel,Xtst);
% acc = mean(double(labeltst == ytst)) * 100;
% fprintf('\nTest Set Accuracy: %f\n', acc);
    
%% Performance parameters for training data

fprintf('\n--------------- Training Data ---------------\n');
letr = length(labeltr);                                % number of training samples
actualtr = zeros(2,letr);                              % actual values of output classes (rows = # of classes, columns = # of training samples) 
predictedtr = zeros(2,letr);                           % predicted values of classes (rows = # of classes, columns = # of training samples)
                                                       % row 1 = will be 1 whenever the sample is of normal class
                                                       % row 2 = will be 1 whenever the sample is of dysphonia class
                                                       % each column is a training sample

for i = 1:2
    actual_index_tr = find(YTrain == i - 1);
    predicted_index_tr = find(labeltr == i - 1);
    actualtr(i,actual_index_tr) = 1;
    predictedtr(i,predicted_index_tr) = 1;
end

performance_stat(actualtr,predictedtr);

% Plotting confusion matrix
figure;
plotconfusion(actualtr,predictedtr);
title('Training Confusion Matrix');

% Plotting Region of Convergence (ROC) curve
figure;
plotroc(actualtr,predictedtr);
title('Training ROC');
 
% Performance parameters for internal cross validating data

fprintf('\n--------------- Cross Validation Data ---------------\n');
letcv = length(labelcv);                                % number of training samples
actualcv = zeros(2,letcv);                              % actual values of output classes (rows = # of classes, columns = # of training samples) 
predictedcv = zeros(2,letcv);                           % predicted values of classes (rows = # of classes, columns = # of training samples)
                                                        % row 1 = will be 1 whenever the sample is of normal class
                                                        % row 2 = will be 1 whenever the sample is of dysphonia class
                                                        % each column is a training sample

for i = 1:2
    actual_index_cv = find(Ycv == i - 1);
    predicted_index_cv = find(labelcv == i - 1);
    actualcv(i,actual_index_cv) = 1;
    predictedcv(i,predicted_index_cv) = 1;
end

performance_stat(actualcv,predictedcv);

% % Plotting confusion matrix
figure;
plotconfusion(actualcv,predictedcv);
title('CV Confusion Matrix');

% % Plotting Region of Convergence (ROC) curve
figure;
plotroc(actualcv,predictedcv);
title('CV ROC');
 
%% Performance parameters for testing data
% 
% fprintf('\n--------------- Cross Validation Data ---------------\n');
% letst = length(labeltst);                                % number of training samples
% actualtst = zeros(2,letst);                              % actual values of output classes (rows = # of classes, columns = # of training samples) 
% predictedtst = zeros(2,letst);                           % predicted values of classes (rows = # of classes, columns = # of training samples)
%                                                         % row 1 = will be 1 whenever the sample is of normal class
%                                                         % row 2 = will be 1 whenever the sample is of dysphonia class
%                                                         % each column is a training sample
% 
% for i = 1:2
%     actual_index_tst = find(ytst == i - 1);
%     predicted_index_tst = find(labeltst == i - 1);
%     actualtst(i,actual_index_tst) = 1;
%     predictedtst(i,predicted_index_tst) = 1;
% end
% 
% performance_stat(actualtst,predictedtst);
% 
% % Plotting confusion matrix
% figure;
% plotconfusion(actualtst,predictedtst);
% title('CV Confusion Matrix');
% 
% % Plotting Region of Convergence (ROC) curve
% figure;
% plotroc(actualtst,predictedtst);
% title('CV ROC');
