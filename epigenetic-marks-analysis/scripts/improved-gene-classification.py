import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.preprocessing import RobustScaler, LabelEncoder
from sklearn.impute import IterativeImputer
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.pipeline import Pipeline

# Constants
DATA_PATH = "../data/features/gene_functionality_features_latest1000all.csv"
FEATURES = [
    "Random_number", "GC_content", "CpG", "GA", "GG", "TA", "TC",
    "lowComplexity_density", "phyloP_max_241w", "phyloP_max_100w",
    "RPKM_tissue", "RPKM_primary.cell", "copy_number", "repeat_free",
    "RNAcode_score", "Max_covariance", "MFE", "RNAalifold_score",
    "Interaction_ave", "gnomAD_SNP_density", "gnomAD_MAF_avg",
    "H3K27ac_AvgSignal", "H3K36me3_AvgSignal", "H3K79me2_AvgSignal",
    "chrm_acc_AvgSignal", "methylome"
]
DATASETS = ['protein-coding-exon2', 'protein-coding-exon3', 'lncrna-exon1', 'lncrna-exon2', 'short-ncrna']

def load_and_preprocess_data(file_path, features, datasets):
    data = pd.read_csv(file_path)
    data = data[data['Dataset'].isin(datasets)]
    X = data[features]
    
    # Simplify dataset mapping
    dataset_mapping = {
        "protein-coding-exon2": "protein-coding",
        "protein-coding-exon3": "protein-coding",
        "lncrna-exon1": "lncrna",
        "lncrna-exon2": "lncrna",
        "short-ncrna": "short-ncrna",
        "protein-exon2-negative-control": "negative-control",
        "protein-exon3-negative-control": "negative-control",
        "lncrna-exon1-negative-control": "negative-control",
        "lncrna-exon2-negative-control": "negative-control",
        "short-ncrna-negative-control": "negative-control"
    }
    data['Dataset'] = data['Dataset'].map(dataset_mapping)
    
    y = data["Dataset"]
    
    return X, y

def create_pipeline():
    return Pipeline([
        ('scaler', RobustScaler()),
        ('imputer', IterativeImputer(max_iter=500)),
    ])

def create_voting_classifier():
    return VotingClassifier(
        estimators=[
            ('lr', LogisticRegression()),
            ('rf', RandomForestClassifier()),
            ('dt', DecisionTreeClassifier(max_depth=10))
        ]
    )

def plot_confusion_matrix(clf, X, y, y_true, name):
    y_pred = cross_val_predict(clf, X, y, cv=5)
    ConfusionMatrixDisplay.from_predictions(
        y_true, y_pred, 
        normalize="true", 
        values_format=".0%", 
        labels=[2, 3, 0, 1],
        display_labels=["protein-coding", "short-ncrna", "lncrna", "negative-control"]
    )
    plt.title(f'{name} Model Confusion Matrix')
    plt.show()

def plot_feature_importance(clf, feature_names, name):
    if hasattr(clf, 'feature_importances_'):
        importances = clf.feature_importances_
        df = pd.DataFrame({'feature': feature_names, 'importance': importances})
        df = df.sort_values('importance', ascending=False)

        plt.figure(figsize=(10, 6))
        plt.barh(df['feature'], df['importance'])
        plt.xlabel('Importance')
        plt.ylabel('Feature')
        plt.title(f'Feature Importance for Gene Type Prediction in {name}')
        plt.show()

def export_decision_tree(clf, feature_names):
    export_graphviz(
        clf,
        out_file="gene_pred_tree.dot", 
        feature_names=feature_names,
        class_names=['lncrna', 'negative-control', 'protein-coding', 'short-ncrna'],
        rounded=True,
        filled=True
    )

def main():
    # Load and preprocess data
    X, y = load_and_preprocess_data(DATA_PATH, FEATURES, DATASETS)
    
    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(y)
    print("Encoded Label equivalence:")
    print(le.transform(["protein-coding", "short-ncrna", "lncrna", "negative-control"]))
    print(le.inverse_transform([0, 1, 2, 3]))
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)
    
    # Create and fit pipeline
    pipeline = create_pipeline()
    X_train_processed = pipeline.fit_transform(X_train)
    X_test_processed = pipeline.transform(X_test)
    
    # Create and train voting classifier
    voting_clf = create_voting_classifier()
    voting_clf.fit(X_train_processed, y_train)
    
    # Evaluate individual classifiers
    for name, clf in voting_clf.named_estimators_.items():
        print(f"{name} accuracy: {clf.score(X_test_processed, y_test):.4f}")
        plot_confusion_matrix(clf, X_test_processed, y_test, y_test, name)
        plot_feature_importance(clf, FEATURES, name)
        
        if name == "Decision Tree Classifier":
            export_decision_tree(clf, FEATURES)
    
    # Evaluate ensemble model
    ensemble_score = voting_clf.score(X_test_processed, y_test)
    print(f"Ensemble Model accuracy: {ensemble_score:.4f}")
    
    # Predict using ensemble model
    prediction = voting_clf.predict(X_test_processed[:1])
    print(f"Prediction for first test sample: {le.inverse_transform(prediction)}")

if __name__ == "__main__":
    main()
