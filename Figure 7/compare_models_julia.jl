using Pkg; Pkg.activate(".")
using FileIO, JLD2
y_data = load("Julia_stochasticity_memRaf1-solns_LHS-1000_tf=15.0_dt=0.1_nreps=5_2021-08-16.jld2") 

using MLJ, CSV, DataFrames
using ModelMiner
X       = CSV.read("LHS_matrix_n=1000_2021-08-16.csv", DataFrame) 
y       = y_data["memraf1_sols"]
results = mine(X,y)

#using PyCall
#@pyimport LazyClassifier
#@pyimport train_test_split
#@pyimport load_breast_cancer

#from lazypredict.Supervised import LazyClassifier
#from sklearn.datasets import load_breast_cancer
#from sklearn.model_selection import train_test_split

#data = load_breast_cancer()
#X = data.data
#y= data.target

#X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=.5,random_state =123)

#clf = LazyClassifier(verbose=0,ignore_warnings=True, custom_metric=None)
#models,predictions = clf.fit(X_train, X_test, y_train, y_test)

#print(models)
using RDatasets
ex_data = dataset("datasets", "iris")
ex_y, ex_X = MLJ.unpack(ex_data, ==(:Species))
ex_results = mine(ex_X, ex_y)