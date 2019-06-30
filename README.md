# Dysphonia-Detection-Using-Gaussian-Mixture-Models
Detecting Dysphonia (voice disorder) using I-vector feature extraction technique and modelling the same using Gaussian Mixture Models.

2 main files are there:

SVMWIvec and SVMIvec

PART 1: SVMWIvec (SVM without using I-Vectors)

Here, navigate to the SVMWIvec folder and click on WIvecGui.fig
A GUI will pop up. We will click all the buttons one by one. 

NOTE: PLEASE NOTE THAT AFTER EVERY BUTTON CLICK, WAIT FOR THE MESSAGE BOX TO SHOW UP

1. Click on 'Read Audio Files' Button and navigate to 'SpeechSamples' folder to select 'TRAIN Dysphonia' and 'TRAIN Normal' folders one by one.
   Wait for a message box to pop up.
2. Click on 'Speech Features' to extract speech features for further processing. Wait till the operation is over.
3. Click on 'Optimization', to train the SVM and find the best hyperparameters for our given data.
4. Click on 'SVM Training' to evaluate our model on SVM Training data. Our training data is broken into Training Set which is used to find best
   hyper parameters, CV set to check which parameters work best and Test Set to evaluate our model on unknown data set.
5. Click on 'SVM New', to check this saved SVM Model and evaluate it's accuracy on our Unknown data set resembling real world data.

NOTE: VERY IMPORTANT: WE HAVE ALREADY FOUND OUT THE BEST MODEL BEFORE HAND. THE ABOVE GUI IS JUST FOR DEMONSTRATION PURPOSE AND ALL THE OPTIMIZATION
AND SPEECH FEATURE EXTRACTION IS JUST FOR UNDERSTANDING PURPOSES. THE 'SVM TRAINING' AND 'SVM NEW' BUTTONS WILL ALWAYS GIVE SAME RESULTS AS WE HAVE LOADED
THE BEST DATA PERMUTATION ALONG WITH THE BEST OPTIMIZED MODEL FOR YOUR UNDERSTANDING. NO REAL TIME PROCESSING IS HAPPENING HERE. 


PART 2: SVMIvec (SVM with I-Vectors)

Here, navigate to the SVMIvec folder and click on IvecGui.fig
A GUI will pop up. We will click all the buttons one by one. 

NOTE: PLEASE NOTE THAT AFTER EVERY BUTTON CLICK, WAIT FOR THE MESSAGE BOX TO SHOW UP

1. Click on 'Read Audio Files' Button and navigate to 'SpeechSamples' folder to select 'TRAIN Dysphonia' and 'TRAIN Normal' folders one by one.
   Wait for a message box to pop up.
2. Click on 'Speech Features' to extract speech features for further processing. Wait till the operation is over.
3. Click 'Calculate I-Vectors' to calculate the 110 I-Vectors using 128 gaussians.
4. Arrange this extracted I-Vector data into a format suitable for SVM Training using 'Arranging Data' button.
5. Click on 'Optimization', to train the SVM and find the best hyperparameters for our given data.
6. Click on 'SVM Training' to evaluate our model on SVM Training data. Our training data is broken into Training Set which is used to find best
   hyper parameters, CV set to check which parameters work best and Test Set to evaluate our model on unknown data set.
7. Click on 'SVM New', to check this saved SVM Model and evaluate it's accuracy on our Unknown data set resembling real world data.

NOTE: VERY IMPORTANT: WE HAVE ALREADY FOUND OUT THE BEST MODEL BEFORE HAND. THE ABOVE GUI IS JUST FOR DEMONSTRATION PURPOSE AND ALL THE OPTIMIZATION
AND SPEECH FEATURE EXTRACTION IS JUST FOR UNDERSTANDING PURPOSES. THE 'SVM TRAINING' AND 'SVM NEW' BUTTONS WILL ALWAYS GIVE SAME RESULTS AS WE HAVE LOADED
THE BEST DATA PERMUTATION ALONG WITH THE BEST OPTIMIZED MODEL FOR YOUR UNDERSTANDING. NO REAL TIME PROCESSING IS HAPPENING HERE. 


Important files to lookout for in both folders:

Folder1: SVMWIvec
1. speech1.m     	: The root file for analyzing speech data without using I-Vectors
2. svmWIvec.m		: SVM training, hyperparameter optimization and final testing on unknown data (sections uncommented as and when required)

Folder1: SVMIvec 
1. speech.m		: The root file for analyzing speech data using I-Vectors
2. i_vector.m		: Contains gmm_em.m, compute_bw_stats.m, train_tv_space.m and extract_ivector.m codes to find I-Vectors. Called from speech.m file.
3. gmm_em.m		: Place where expectation and maximization done on the incoming data and generates a UBM model (returns UBM)
4. compute_bw_stats.m	: Here, mapping of each sample to UBM model takes place. It generates supervector using 0th and 1st order statistics. 
                          It uses, Baum Welsh Algorithm hence the name BW_stats.
5. train_tv_space.m	: Here we find the total variability in the model.
6. extract_ivector.m	: It is used to extract I-Vectors using the above TV subspace.
7. data_arrange.m	: After I-Vector extraction, it is necessary for the data to be arranged in format useful for training and predictions. Contains wccn.m
8. wccn.m		: We need to apply Cholesky whitening after I-Vector extraction to remove within session variability.
9. svmIvec.m		: SVM training, hyperparameter optimization and final testing on unknown data (sections uncommented as and when required)

Some common codes in both the folders which are required.
1. results.m		: To compute prediction results and find accuracies and plot data whenever relevant.
2. performance_stat.m	: To find performance paramters like TPR, TNR, FAR, FRR, F1 statistic, etc.
3. performance.m	: To plot ROC and confusion matrices
4. featureNormalize.m	: Normalizing the features
5. randomize.m		: To randomnly permutate the data matrices as and when required and break them into training, CV and test sets in 70:20:10 ratio.
