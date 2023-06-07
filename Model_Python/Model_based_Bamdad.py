import os
import numpy as np
import scipy.io
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GridSearchCV
from geneticalgorithm import geneticalgorithm as ga
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import h5py

# Define the transfer function
def transfer_function(x, a, b):
    return a - np.exp(-x/b)**2

# Define function to load MATLAB file with nested structs
def load_matlab_struct(file_path, struct_path):
    """
    Args:
        file_path (str): The path to the .mat file.
        struct_path (str): The path to the struct in the .mat file, using dot notation to specify nested structs.

    Returns:
        A dictionary containing the struct fields and values.
    """
    # Load the data from the .mat file
    with h5py.File(file_path, 'r') as file:
        struct_data = file[struct_path]
    # Convert the nested structs to dictionaries
        if isinstance(struct_data, h5py.Group):
            struct_data = _convert_h5py_group(struct_data)

    return struct_data

def _convert_h5py_group(group):
    """
    Recursively converts an h5py Group object to a dictionary.
    """
    dict_ = {}
    for key, val in group.items():
        if isinstance(val, h5py.Group):
            dict_[key] = _convert_h5py_group(val)
        else:
            dict_[key] = val[()]
    return dict_import glob

from scipy.signal import iirnotch, filtfilt
# Create empty lists to store the features and labels
features = []
labels = []
amplitudes = []
baselines = []
gradients = []
extension = '.bin'
file_pattern = '_QUAT.mat'

# Define the time points to extract data from
baseline_time = 2.0
gradient_start_time = 4.0
gradient_end_time = 6.0
output_time = 9.0

# Loop over files in folder
folder_path = "/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/single_subject_dataset_10%"
for subfolder_name in os.listdir(folder_path):
    subfolder_path = os.path.join(folder_path, subfolder_name)
    # Ask the optimal channel
    channel = None
    while channel is None:
        print('Now is processing: ', subfolder_path)
        channel_str = input("Enter the optimal channel: ")
        try:
            channel = int(channel_str)
            if channel < 1 or channel > 16:
                raise ValueError()
        except ValueError:
            print("Invalid channel. Please enter a number between 1 and 16.")
            channel = None

    channel -= 1

    if os.path.isdir(subfolder_path):
        # Get a list of sorted file names
        #sorted_files = sorted([f for f in os.listdir(subfolder_path) if f.endswith(extension)])
        to_scan = subfolder_path + '/*_QUAT*.mat'
        sorted_files = sorted(glob.glob(to_scan))
        # Load the additional data from the *.mat file
        with h5py.File(os.path.join(subfolder_path, sorted_files[0][:-3] + 'mat'), 'r') as f:
        # Get the trial stimulation amplitudes for all files
            #for key in f.keys():
                #print(key)
            vib_amps = f['RecInfo']['Experiment']['Order'][:]

        # Loop over the sorted file names and load the raw EMG from the BIN file
        for nfile, file_name in enumerate(sorted_files):
            # Input 1 - stimulation amplitude for the current file
            with h5py.File(os.path.join(subfolder_path, file_name[:-3] + 'mat'), 'r') as f:
                if nfile == 1:
                    vib_amp = vib_amps[0]
                else:
                    vib_amp = vib_amps[nfile-1]
                # Get the TimeStamp data
                Time_Stamps = f['RecInfo']['Data']['TimeStamp'][:]
                # Remove the time offset
                Time = Time_Stamps - np.ones(len(Time_Stamps)) * Time_Stamps[0, 0]

                # Get EMVC
                EMVC = f['RecInfo']['Calibration']['EMVC'][:]

            # Load the raw EMG from the BIN file
            bin_name = file_name[:-3] + 'bin'
            with open(os.path.join(subfolder_path, bin_name), 'rb') as f:
                raw_emg = np.fromfile(f, dtype=np.int16)
                raw_emg = raw_emg.reshape(-1, 16).T

                # Filter noise
                qfactor = 10
                fe = 116  # Assign peak frequency
                n_harmonics = 2
                fs = 1000  # system at 1k hz

                if fe > 0:
                    frequencies = np.arange(1, n_harmonics + 1) * fe
                    for freq_to_notch in frequencies:
                        wo = freq_to_notch / (fs / 2)
                        bw = wo / qfactor
                        b, a = iirnotch(wo, bw)
                        raw_emg = filtfilt(b, a, raw_emg)

                # Filter noise *2
                qfactor = 10
                fe = 49  # Assign peak frequency
                n_harmonics = 2
                fs = 1000  # system at 1k hz

                if fe > 0:
                    frequencies = np.arange(1, n_harmonics + 1) * fe
                    for freq_to_notch in frequencies:
                        wo = freq_to_notch / (fs / 2)
                        bw = wo / qfactor
                        b, a = iirnotch(wo, bw)
                        raw_emg = filtfilt(b, a, raw_emg)

                # EMG raw data at the optimal channel
                rmsWindow = 500
                BufSize = 40
                nSamples = raw_emg.shape[1]
                nWindows = int((nSamples / BufSize) - (rmsWindow / BufSize)) #number of windows
                samplect = 0
                nChannels = 16
                rmsEMG = np.zeros(nWindows) #number of channels - 16 - rmsEMG array size

                for n in range(1, nWindows + 1):
                    if n == 1:
                        DataI = raw_emg[channel, :rmsWindow]
                    else:
                        start_idx = (n - 1) * BufSize + 1
                        end_idx = start_idx + rmsWindow - 1
                        DataI = raw_emg[channel, start_idx:end_idx]
                    rmsEMG[samplect] = np.sqrt(np.mean(np.square(DataI)))
                    samplect += 1

                # Normalise data
                Data = np.zeros(nWindows)
                for i in range(nWindows):
                    Data [i] = (rmsEMG[i] - EMVC[channel, 1]) / (EMVC[channel, 0] - EMVC[channel, 1])

                # Input 2 - baseline data (normalized EMG data at 2s)
                Time_shape = Time.shape
                Data_shape = Data.shape
                length = min(Time_shape[1], Data_shape[0])
                Time_slice = Time[:, :length]
                Data_slice = Data[:length]
                Time_bcast = Time_slice.reshape(length, 1)
                Data_bcast = Data_slice.reshape(-1, 1)
                pairs = np.hstack((Time_bcast, Data_bcast))
                index = np.abs(pairs[:, 0] - 2).argmin()
                baseline_data = pairs[index, 1]

                # Input 3 - gradient between 4s and 6s - use delta EMG instead as more accurate
                index_after_stimulation = np.abs(pairs[:, 0] - 6).argmin()
                after_stimulation_data = pairs[index_after_stimulation, 1]
                """
                extracted_data = pairs[(pairs[:, 0] >= 4) & (pairs[:, 0] <= 6)]
                emg_data_for_gradient = extracted_data[:, 1]
                time_values_for_gradient = extracted_data[:, 0]
                reg = LinearRegression().fit(time_values_for_gradient.reshape(-1, 1), emg_data_for_gradient)
                """
                gradient = (after_stimulation_data - baseline_data)/1.5

                # Append the input features to the corresponding lists
                features.append([vib_amp, gradient, baseline_data])
                amplitudes.append(vib_amp)
                baselines.append(baseline_data)
                gradients.append(gradient)

                # Fit the transfer function to the data and append the output
                # Define the threshold for the covariance matrix
                cov_threshold = 0.5 #tunable

                extracted_plateau = pairs[(pairs[:, 0] >= 6) & (pairs[:, 0] <= 9)]
                emg_data_for_transfer = extracted_plateau[:, 1]
                time_values_for_transfer = extracted_plateau[:, 0]
                popt, pcov = curve_fit(transfer_function, time_values_for_transfer, emg_data_for_transfer)
                """
                if np.all(np.diag(pcov) > cov_threshold):
                    # peak activity output at 9s
                    a_opt = popt[0]
                    b_opt = popt[1]
                    peak_value = transfer_function(9, a_opt, b_opt)
                else:
                    continue
                """
                a_opt = popt[0]
                b_opt = popt[1]
                peak_value = transfer_function(9, a_opt, b_opt)
                labels.append(peak_value)
# Convert the lists to numpy arrays
X = np.array(features, dtype=object)
y = np.array(labels)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

print(X_train.shape)
print(y_train.shape)

import csv
import numpy as np

# Define the file name and path
file_name = 'data.csv'
file_path = '/Users/xinyaoyo/MATLAB/Projects/UROP/Codes/Data_Jupyter'

# Combine the feature set and label set into one array
data = np.column_stack((X_train, y_train))
header = ['feature1', 'feature2', 'feature3', 'label']

# Save the data to a CSV file
with open(file_path + '/' + file_name, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    for row in data:
        row_str = [str(r) for r in row]
        writer.writerow(row_str)

import numpy as np
from sklearn.model_selection import KFold, cross_val_score
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from tpot import TPOTClassifier

# Make sure Structure
X_train = np.array([np.concatenate(x) for x in X_train])
X_train = X_train.astype(np.float)

y_train = np.array(y_train)
y_train = y_train.astype(np.float)

# Define the classifiers and their parameter spaces for optimization with GA
svm_params = {
    'C': np.logspace(-4, 4, 10),
    'kernel': ['linear', 'poly', 'rbf', 'sigmoid'],
    'degree': [2, 3, 4],
    'gamma': ['scale', 'auto'] + list(np.logspace(-4, 4, 10)),
    'coef0': np.linspace(-1, 1, 21),
}

rf_params = {
    'n_estimators': [10, 50, 100, 500],
    'criterion': ['gini', 'entropy'],
    'max_depth': [None] + list(range(5, 21)),
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['sqrt', 'log2', None],
}

# Create the classifiers
svm = SVC()
rf = RandomForestClassifier()

# Create the TPOT optimizer
tpot = TPOTClassifier(generations=5, population_size=20, cv=5, random_state=42, verbosity=2)

# Create the K-Fold cross-validator
kfold = KFold(n_splits=10, shuffle=True, random_state=42)

# Evaluate each classifier using K-Fold cross-validation
svm_scores = cross_val_score(svm, X_train, y_train, cv=kfold)
rf_scores = cross_val_score(rf, X_train, y_train, cv=kfold)
tpot_scores = cross_val_score(tpot, X_train, y_train, cv=kfold)

# Print the mean accuracy for each classifier
print(f'SVM mean accuracy: {np.mean(svm_scores):.3f}')
print(f'RF mean accuracy: {np.mean(rf_scores):.3f}')
print(f'TPOT mean accuracy: {np.mean(tpot_scores):.3f}')

# Optimize the hyperparameters of the classifiers using TPOT
svm_params_opt = tpot.fit(X_train, y_train).get_params()['pipeline'].get_params()
rf_params_opt = tpot.fit(X_train, y_train).get_params()['pipeline'].get_params()
# Print the hyperparameters of the optimized classifiers
print(f'Optimized SVM parameters: {svm_params_opt}')
print(f'Optimized RF parameters: {rf_params_opt}')

# Train the optimized classifiers and test their accuracy on the test set
X_test = np.array([np.concatenate(x) for x in X_test])
X_test = X_test.astype(np.float)
y_test = np.array(y_test)
y_test = y_test.astype(np.float)

svm_opt = SVC(**svm_params_opt)
svm_opt.fit(X_train, y_train)
svm_opt_score = svm_opt.score(X_test, y_test)

rf_opt = RandomForestClassifier(**rf_params_opt)
rf_opt.fit(X_train, y_train)
rf_opt_score = rf_opt.score(X_test, y_test)

# Print the accuracy of the optimized classifiers on the test set
print(f'Optimized SVM accuracy on test set: {svm_opt_score}')
print(f'Optimized RF accuracy on test set: {rf_opt_score}')

from numpy import float64
# Method 2 - blood test - link: https://github.com/ritabratamaiti/Blooddonorprediction/blob/master/script.py
# Compare SVM, Perceptron, k-NN, Naive Bayes, and decision trees TPOTClassifier classification methods
from sklearn.svm import SVC
from sklearn.linear_model import Perceptron
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.model_selection import KFold
from sklearn import tree
import numpy as np
from tpot import TPOTClassifier

X_train = X_train.astype(np.float64)
y_train = y_train.astype(np.float64)
tpot = TPOTClassifier(generations=20, population_size=50, verbosity=2)#large generations for small datasets, larger population means better results
tpot.fit(X_train, y_train)

#classifiers
clf_tree = tree.DecisionTreeClassifier()
clf_svm = SVC()
clf_perceptron = Perceptron()
clf_KNN = KNeighborsClassifier()
clf_nb = BernoulliNB()

list_clf = [clf_tree, clf_svm, clf_perceptron, clf_KNN, clf_nb, tpot.fitted_pipeline_]

#K-fold cross-validation
kf = KFold(n_splits=5)
kf.get_n_splits(X)
c = 1
for clfs in list_clf:
    print(c)
    c += 1
    a = 0
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        clfs.fit(X_test, y_test)
        print(clfs.score(X_train, y_train))
        a += clfs.score(X_train, y_train)
    a = a/5
    print("Average=",a,"\n")
    print(clfs,"\n")