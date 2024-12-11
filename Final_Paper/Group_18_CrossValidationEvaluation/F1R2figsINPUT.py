import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.datasets import load_diabetes
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import PredictionErrorDisplay
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE

# !!!
# Here, we get the R^2 (coefficient of determination) values for training and testing data.
# "R-squared is a statistical measure that indicates how much of the variation of a dependent variable is explained by an independent variable in a regression model."
# In essence, the higher the score, the better predictive performance. 
# !!!


filepath = input("Enter path to the CSV: ")
data = pd.read_csv(filepath)
print(" ")
print(" ")
print("R^2 Evaluation:")
print(" ")

X = data.iloc[:, :-1].values  # All columns except the last one
y = data.iloc[:, -1].values   # The last column


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

model = LinearRegression().fit(X_train, y_train)
print("Train-Test Split, Linear Regression R^2 Training Data: ")
trainingScore = model.score(X_train, y_train)
print(trainingScore)
print("Train-Test Split, Linear Regression R^2 Testing Data: ")
testingScore = model.score(X_test, y_test)
print(testingScore)

# PLOT #
y_pred = model.predict(X_train)
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
PredictionErrorDisplay.from_predictions(
    y[:len(y_pred)],
    y_pred=y_pred,
    kind="actual_vs_predicted",
    subsample=100,
    ax=axs[0],
    random_state=0,
)
axs[0].set_title("Actual vs. Predicted values")
PredictionErrorDisplay.from_predictions(
    y[:len(y_pred)],
    y_pred=y_pred,
    kind="residual_vs_predicted",
    subsample=100,
    ax=axs[1],
    random_state=0,
)
axs[1].set_title("Residuals vs. Predicted Values")
fig.suptitle("Plotting Train-Test Split predictions")
plt.tight_layout()
plt.show()


print(" ")
print(" ")
# K-Folds (5) #

k = 5
kf = KFold(n_splits=k, shuffle=True, random_state=42)
model = LinearRegression()
scores = cross_val_score(model, X, y, cv=kf, scoring='r2')

average_r2 = np.mean(scores) 
print("K-folds, Linear Regression: ")
print(f"R² Score for each fold: {[round(score, 4) for score in scores]}")
print(f"Average R² across {k} folds: {average_r2:.2f}")

# PLOT #
y_pred = cross_val_predict(model, X, y, cv=kf)
fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
PredictionErrorDisplay.from_predictions(
    y,
    y_pred=y_pred,
    kind="actual_vs_predicted",
    subsample=100,
    ax=axs[0],
    random_state=0,
)
axs[0].set_title("Actual vs. Predicted values")
PredictionErrorDisplay.from_predictions(
    y,
    y_pred=y_pred,
    kind="residual_vs_predicted",
    subsample=100,
    ax=axs[1],
    random_state=0,
)
axs[1].set_title("Residuals vs. Predicted Values")
fig.suptitle("Plotting cross-validated predictions")
plt.tight_layout()
plt.show()

print(" ")
print(" ")
print("F1 Evaluation:")
print(" ")
# best threshold using percentiles
best_threshold = None
best_f1 = 0
thresholds = np.percentile(y, range(5, 95, 5))  #
#percentiles from 5% to 95%

for t in thresholds:
    yclass = (y > t).astype(int)
    f1_scores = []
    for train_index, test_index in kf.split(X, yclass):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = yclass[train_index], yclass[test_index]
        model = LogisticRegression()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        f1_scores.append(f1_score(y_test, y_pred))
    avg_f1 = np.mean(f1_scores)
    if avg_f1 > best_f1:
        best_f1 = avg_f1
        best_threshold = t

print(f"Best threshold: {best_threshold}, Best F1 Score: {best_f1:.2f}")
yclass = (y > best_threshold).astype(int)

scaler = StandardScaler()
X = scaler.fit_transform(X)
smote = SMOTE(random_state=42)
Xresample, yresample = smote.fit_resample(X, yclass)


# Train Test Split F1
X_train, X_test, y_train, y_test = train_test_split(Xresample, yresample, test_size=0.2, random_state=0)
model = LogisticRegression().fit(X_train, y_train)

ytrainpred = model.predict(X_train)
ytestpred = model.predict(X_test)
trainf1 = f1_score(y_train, ytrainpred)
testf1 = f1_score(y_test, ytestpred)

print("Train-Test Split, Logistic Regression (F1):")
print(f"F1 Score (Training Data): {trainf1:.4f}")
print(f"F1 Score (Testing Data): {testf1:.4f}")


print(" ")
print(" ")


# K-Fold Cross Validation F1 (5)
k = 5
kf = KFold(n_splits=k, shuffle=True, random_state=42)
model = LogisticRegression()

f1_scores = []

for train_index, test_index in kf.split(X):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = yclass[train_index], yclass[test_index]
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    f1_scores.append(float(f1_score(y_test, y_pred)))  # Convert to Python float

average_f1 = np.mean(f1_scores)
print("\nK-folds, Logistic Regression (F1, Representing Testing Sets):")
print(f"F1 Score for each fold: {[round(score, 4) for score in f1_scores]}")
print(f"Average F1 Score across {k} folds: {average_f1:.4f}")
