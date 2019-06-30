function data_arrange(devIVs,~)

[ivec,class,samples] = size(devIVs);

%% Initializing necessary parameters

t = zeros(ivec,class * samples);             
c1 = zeros(ivec,samples);                   % I-Vectors for class 1
c2 = zeros(ivec,samples);                   % I-sVectors for class 2

%% Extracting data of all classes into a vector

parfor i = 1:samples * class
    t(:,i) = devIVs(:,i);
end

%% Seperating data according to individual classes and creating targets

for j = 1:samples
    c1(:,j) = t(:,2 * j - 1);
    c2(:,j) = t(:,2 * j);
end

iVector = [c1'; c2'];
y = ones(samples * class,1);
y(samples + 1:end,:) = 0;

%% Randomly permutate the extracted I-Vectors

sel = randperm(size(iVector, 1));          % holds indices of randomly arranged data values
iVector = iVector(sel,:);                  % arrange data depending upon above indices
y = y(sel,:);                              % arrange targets depending upon above indices

%% Performing Cholesky whitening on the extracted I-Vectors

alpha = 0.1;                               % Cholesky whitening co-efficient
[X,~] = wccn(iVector,y,alpha);

%% Organizing data for further processing

% setting percentages in which data has to be broken
[h,~] = size(X);
p1 = ceil(0.7 * h);
p2 = p1 + ceil(0.2 * h);

% creating file handlers
filename = 'datat.mat';
fl = matfile(filename,'Writable',true);

% training data
fl.Xtr = X(1:p1,:);
fl.ytr = y(1:p1,:);

% cross validation data
fl.Xcv = X(p1 + 1:p2,:);
fl.ycv = y(p1 + 1:p2,:);

% test data
fl.Xtst = X(p2 + 1:end,:);
fl.ytst = y(p2 + 1:end,:);

% total data
fl.X = X;
fl.y = y;

%% Saving data

save indexest.mat sel;
end