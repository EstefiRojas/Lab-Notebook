from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.tree import export_graphviz
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import ConfusionMatrixDisplay

import graphviz

import pandas as pd

import numpy as np

import matplotlib.pyplot as plt


# Variables
select_datasets = ['protein-coding-exon2','protein-coding-exon3','protein-exon2-negative-control','protein-exon3-negative-control']
select_datasets = ['lncrna-exon1','lncrna-exon2','lncrna-exon1-negative-control','lcnrna-exon2-negative-control']
select_datasets = ['short-ncrna','short-ncrna','short-ncrna-negative-control','short-ncrna-negative-control']
select_datasets = ['protein-coding-exon2','protein-coding-exon3','lncrna-exon1','lncrna-exon2','short-ncrna']


# 1. Load features
data = pd.read_csv("../data/features/gene_functionality_features_latest1000all.csv", sep=",")

# Features of interest
#select_features = ["Random number",
#                   "GC content",
#                   "CpG","GA","GG","TA","TC",
#                   "low_complexity_density","phyloP max_241w","phyloP max_100w",
#                   "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max_covariance",
#                   "MFE","RNAalifold_score","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF_avg",
#                   "H3K27ac_AvgSignal","H3K36me3_AvgSignal","H3K79me2_AvgSignal","chrm_acc_AvgSignal","methylome"]
# Features of interest
select_features = ["Random number",
                   "GC content",
                   "CpG","GA","GG","TA","TC",
                   "low_complexity_density","phyloP max_241w","phyloP max_100w",
                   "RPKM_tissue","RPKM_primary cell","Copy number","Repeat free","RNAcode","Max covariance",
                   "MFE","RNAalifold","Interaction_ave","gnomAD_SNP_density","gnomAD_MAF",
                   "H3K27ac","H3K36me3","H3K79me2","chromatin_acc","methylome"]

# 1.1 Preprocess data
data = data.dropna() # Drop rows with missing values
#data = data[((data['Dataset'] == select_datasets[0]) | 
#             (data['Dataset'] == select_datasets[1]) | 
#             (data['Dataset'] == select_datasets[2]) | 
#             (data['Dataset'] == select_datasets[3]))]


# 2. Features (X) and target variable (y)
X = data.drop(["Functional","Dataset"], axis=1)
X = X[select_features]

y = data["Functional"]
#data.loc[data['Dataset'] == "protein-coding-exon2", 'Dataset'] = "protein-coding"
#data.loc[data['Dataset'] == "protein-coding-exon3", 'Dataset'] = "protein-coding"
#data.loc[data['Dataset'] == "protein-exon2-negative-control", 'Dataset'] = "negative-control"
#data.loc[data['Dataset'] == "protein-exon3-negative-control", 'Dataset'] = "negative-control"

#data.loc[data['Dataset'] == "lncrna-exon1", 'Dataset'] = "lncrna"
#data.loc[data['Dataset'] == "lncrna-exon2", 'Dataset'] = "lncrna"
#data.loc[data['Dataset'] == "lncrna-exon1-negative-control", 'Dataset'] = "negative-control"
#data.loc[data['Dataset'] == "lncrna-exon2-negative-control", 'Dataset'] = "negative-control"

#data.loc[data['Dataset'] == "short-ncrna-negative-control", 'Dataset'] = "negative-control"

#y = data["Dataset"]

# Preprocess X
transformer = RobustScaler().fit(X)
X = transformer.transform(X)

iterative_imputer = IterativeImputer(max_iter=500)
X = iterative_imputer.fit_transform(X)
#X = round(pd.DataFrame(iterative_imputed_data, columns=X.columns), 3)

# Preprocess y
le = preprocessing.LabelEncoder()
#le.fit(["protein-coding", "short-ncrna", "lncrna", "negative-control"])
le.fit(["0", "1"])
print("Encoded Label equivalence:")
#print(le.transform(["protein-coding", "short-ncrna", "lncrna", "negative-control"]))
print(le.transform(["0", "1"]))
#print(le.inverse_transform([0, 1, 2, 3]))
print(le.inverse_transform([0, 1]))
y = le.transform(y)


# 3. Split data into training and testing
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)


# 4. Define a hard voting classifier
voting_clf = VotingClassifier(
	estimators=[
		('Linear Regression Classifier', LogisticRegression()),
		('Random Forest Classifier', RandomForestClassifier()),
		('Decision Tree Classifier', DecisionTreeClassifier(max_depth=10))
	]
)


# 5. Train the models
voting_clf.fit(X_train, y_train)


# 6. Print precision scores
for name, clf in voting_clf.named_estimators_.items():
	print(name, "=", clf.score(X_test, y_test))
	#print(name, "=", cross_val_score(clf, X_test, y_test, cv = 3, scoring = "accuracy"))
	y_test_pred = cross_val_predict(clf, X_test, y_test, 
									cv = 5)
	print(y_test_pred)
	ConfusionMatrixDisplay.from_predictions(y_test, y_test_pred, 
											normalize = "true", 
											values_format = ".0%",
											labels = [0,1], 
											#labels = [2,3,0,1],
											display_labels=["No", "Yes"],
											#display_labels=["protein-coding", "short-ncrna", "lncrna", "negative-control"],
											)
	plt.title(name + ' Model Confusion Matrix')
	plt.show()
	if(name=="Decision Tree Classifier" or name=="Random Forest Classifier"):
		feature_names = select_features
		importances = clf.feature_importances_

		if(name=="Decision Tree Classifier"):
			dot_data = export_graphviz(
				clf,
				out_file="gene_pred_tree.dot", 
				feature_names=feature_names,
				#class_names=['lncrna', 'negative-control', 'protein-coding', 'short-ncrna'],
				class_names=['No', 'Yes'],
				rounded=True,
				filled=True)
			#graph = graphviz.Source(dot_data) 
			#graph.render("gene_pred_tree")

		df = pd.DataFrame({'feature': feature_names, 'importance': importances})
		df = df.sort_values('importance', ascending=False)

		plt.figure(figsize=(10, 6))
		plt.barh(df['feature'], df['importance'])
		plt.xlabel('Importance')
		plt.ylabel('Feature')
		plt.title('Feature Importance for Function Prediction in ' + name)
		plt.show()


# 7. Predict test set using ensemble model
voting_clf.predict(X_test[:1])


# 8. Get score for ensemble model
print("Hard Voting ", voting_clf.score(X_test, y_test))

