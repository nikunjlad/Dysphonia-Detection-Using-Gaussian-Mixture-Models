function [acc,TP,FP,precision,TPR,TNR,FAR,FRR,Fmeasure,AUC] = results(label,y,setname,flag)

acc = mean(double(label == y)) * 100;
fprintf('\n%s Set Accuracy: %f\n', setname, acc);
le = length(label);                                % number of training samples
actual = zeros(2,le);                              % actual values of output classes (rows = # of classes, columns = # of training samples) 
predicted = zeros(2,le); 

for i = 1:2
    actual_index = y == i - 1;
    predicted_index = label == i - 1;
    actual(i,actual_index) = 1;
    predicted(i,predicted_index) = 1;
end

%% Calculating performance parameters

[TP,FP,precision,TPR,TNR,FAR,FRR,Fmeasure] = performance_stat(actual,predicted);

%% Finding AUC 
if flag == 1
    [~,~,~,AUC,~] = perfcurve(y,label,1);
    fprintf('\n AUC: %f\n', AUC);
else
    AUC = 0;
end

%% Plotting confusion matrix
figure;
plotconfusion(actual,predicted);
title(strcat(setname,' Confusion Matrix'));

%% Plotting Region of Convergence (ROC) curve
figure;
plotroc(actual,predicted);
title(strcat(setname,' ROC'));

end