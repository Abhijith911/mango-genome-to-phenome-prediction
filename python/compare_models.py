import pandas as pd
import joblib
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from sklearn.svm import SVR
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import json

# -----------------------------
# SNP ENCODING
# -----------------------------
def encode_snp_column(series, snp_name):

    parts = snp_name.split("_")
    if len(parts) < 5:
        return series

    ref = parts[3]
    alt = parts[4]

    het_codes = {"R", "Y", "S", "W", "K", "M"}

    def convert(x):
        if pd.isna(x):
            return np.nan

        x = str(x)

        if x == ref:
            return 0
        elif x == alt:
            return 2
        elif x in het_codes:
            return 1
        elif x == "N":
            return np.nan

        return np.nan

    return series.apply(convert)


def encode_dataset(df):

    df = df.copy()

    snp_cols = [
        c for c in df.columns
        if c not in ["Accession", "SeedThickness"]
    ]

    for col in snp_cols:
        df[col] = encode_snp_column(df[col], col)

    return df


# -----------------------------
# DATASETS
# -----------------------------
datasets = ["Top100", "Top500", "Top1000"]

cv = KFold(n_splits=5, shuffle=True, random_state=42)

results = []

# -----------------------------
# LOOP
# -----------------------------
for ds in datasets:

    print(f"\nRunning {ds}")

    train_df = pd.read_csv(f"SeedThickness_{ds}_Train.csv")
    test_df = pd.read_csv(f"SeedThickness_{ds}_Test.csv")

    train_df = encode_dataset(train_df)
    test_df = encode_dataset(test_df)

    X_train = train_df.drop(columns=["Accession", "SeedThickness"])
    X_test = test_df.drop(columns=["Accession", "SeedThickness"])

    y_train = train_df["SeedThickness"]
    y_test = test_df["SeedThickness"]

    # -------------------------
    # IMPUTE
    # -------------------------
    imputer = SimpleImputer(strategy="most_frequent")

    X_train = pd.DataFrame(imputer.fit_transform(X_train), columns=X_train.columns)
    X_test = pd.DataFrame(imputer.transform(X_test), columns=X_test.columns)
    joblib.dump(imputer, f"imputer_{ds}.pkl")

    # -------------------------
    # MODELS (ONLY 2 YOU NEED)
    # -------------------------
    models = {

        "LightGBM": LGBMRegressor(
            n_estimators=500,
            learning_rate=0.05,
            random_state=42
        ),

        "ElasticNet": Pipeline([
            ("scaler", StandardScaler()),
            ("model", ElasticNet(
                alpha=0.1,
                l1_ratio=0.5,
                max_iter=10000
            ))
        ])
    }

    # -------------------------
    # TRAIN + EVAL
    # -------------------------
    for model_name, model in models.items():

        cv_scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="r2")

        model.fit(X_train, y_train)
        joblib.dump(model, f"{model_name}_{ds}.pkl")

        train_pred = model.predict(X_train)
        test_pred = model.predict(X_test)

        results.append({
            "Dataset": ds,
            "Model": model_name,
            "Train_R2": r2_score(y_train, train_pred),
            "CV_R2_Mean": np.mean(cv_scores),
            "Test_R2": r2_score(y_test, test_pred),
            "RMSE": np.sqrt(mean_squared_error(y_test, test_pred)),
            "MAE": mean_absolute_error(y_test, test_pred)
        })

        print(ds, model_name, round(r2_score(y_test, test_pred), 4))

# -----------------------------
# SAVE RESULTS
# -----------------------------
results_df = pd.DataFrame(results)
results_df = results_df.sort_values("Test_R2", ascending=False)

results_df.to_excel("SeedThickness_Model_Comparison.xlsx", index=False)

print("\nDONE")

for ds in ["Top100", "Top500", "Top1000"]:
    train_df = pd.read_csv(f"SeedThickness_{ds}_Train.csv")
    features = [c for c in train_df.columns if c not in ["Accession", "SeedThickness"]]
    with open(f"features_{ds}.json", "w") as f:
        json.dump(features, f)
print("Feature JSONs saved.")