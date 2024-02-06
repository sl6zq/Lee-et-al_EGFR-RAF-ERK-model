import h5py

f = h5py.File("Julia_stochasticity_memRaf1-solns_LHS-1000_tf=15.0_dt=0.1_nreps=5_2021-08-16.jld2", "r")
f["stoch_sols"].value


from lazypredict.Supervised import LazyClassifier
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split

data = load_breast_cancer()
X = data.data
y= data.target

X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=.5,random_state =123)

clf = LazyClassifier(verbose=0,ignore_warnings=True, custom_metric=None)
models,predictions = clf.fit(X_train, X_test, y_train, y_test)

print(models)