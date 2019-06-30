
id = acc(:,1) == max(acc(:,1));
trval = max(acc(:,1));
tstval = acc(id,2);
tstval = min(tstval);
val = [trval,tstval];
id = find(acc(:,1) == trval & acc(:,2) == tstval, 1, 'last' );

[m1,n1] = size(a);
for i = 1:m1
    for j = 1:n1
        if isequaln(a{i,j},val)
            lasso_component = i * 2 - 1;
            no_of_ivectors = 110 + (j * 5) - 5;
        end
    end
end

figure;
plot(acc(:,1));
hold on;
stem(id,ones(length(id),1) * trval);
legend('Training accuracies','Optimal parameters');
xlabel('No of probable samples');
ylabel('Training Accuracies');
text(id + 5,trval,['Lasso: ' num2str(lasso_component),', I-Vectors: ' num2str(no_of_ivectors)]);
text(id + 5,trval / 2,['Training: ' num2str(trval) '%, Testing: ' num2str(tstval) '%']);

% save C:\Users\etrx\Desktop\shipnik\Dysphonia\gauss\Training_using_Optimal_parameters\optimum_params.mat lasso_component no_of_ivectors;

% save G:\shipnik\Dysphonia\gauss\Training_using_Optimal_parameters\optimum_params.mat lasso_component no_of_ivectors;