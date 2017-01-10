import sys
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import load_svmlight_file
from numpy import around

argvs = sys.argv
argc  = len(argvs)

model_path            = argvs[1]
test_feature_vec_path = argvs[2]

with open(model_path, mode='rb') as f:
    clf = pickle.load(f)

# vector laoding
X, y = load_svmlight_file(test_feature_vec_path)

# prediction
output_proba = clf.predict_proba(X)

# output sequentially
for i in output_proba:
    prob = around(i[1], 2)
    print prob
    
