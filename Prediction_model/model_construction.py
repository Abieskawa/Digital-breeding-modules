#!/usr/bin/env python3
"""
Step 3: Model Construction

Inputs:
  - Model-ready train/test CSVs produced by Step 1 (model_prep):
      output_dir/Model_result/{file_fold_index}_{pc}PCs_{n}SNPs/
        train_selected_snps_*.csv
        test_selected_snps_*.csv

Outputs:
  - Metrics CSVs (per model/config)
  - Saved models (NEW) under output_dir/Model_result/Saved_models/

Cross-validation orchestration is intentionally NOT implemented here
(you said CV stays in the main run_prediction_pipeline.py).
"""
import os
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np
import pandas as pd
from Utils.utils import _resolve_outdir, parse_int_list

from sklearn import metrics, ensemble, svm
from sklearn.metrics import confusion_matrix
from xgboost import XGBClassifier

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.init as init
from torch.utils.data import DataLoader, TensorDataset


try:
    from sklearnex import patch_sklearn, unpatch_sklearn
except Exception:
    patch_sklearn = None
    unpatch_sklearn = None

try:
    import joblib
except Exception:
    joblib = None


class nnModel(nn.Module):
    def __init__(self, num_i, num_h, num_o):
        super(nnModel, self).__init__()
        self.fc1 = nn.Linear(num_i, num_h)
        self.fc2 = nn.Linear(num_h, num_h)
        self.fc3 = nn.Linear(num_h, num_o)

        init.xavier_normal_(self.fc1.weight)
        init.xavier_normal_(self.fc2.weight)
        init.xavier_normal_(self.fc3.weight)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x


class ModelConstruction:
    def __init__(self, config: dict):
        self.config = config

        self.species_name = self.config.get('species_name')
        self.phenotype_column = self.config.get('phenotype_column', 'Phenotype')
        self.data_type = self.config.get('data_type', 'binary')

        self.output_dir = _resolve_outdir(
            self.config,
            key="prediction_output_dir",
            default="results",
            ensure_dir=True,
        )

        self.model_output_dir = _resolve_outdir(
            base_outdir=self.output_dir,
            subdir="Model_result",
            ensure_dir=True,
        )

        # Defaults from config
        self.models = [m.strip() for m in self.config.get('models', 'SVM,RandomForest,XGBoost').split(',') if m.strip()]
        self.top_n_snps_list = parse_int_list(self.config.get('top_n_snps_list', '10'))
        self.gwas_n_pcs_list = parse_int_list(self.config.get('gwas_n_pcs_list', '')) or \
                              parse_int_list(self.config.get('gwas_n_pcs', '')) or \
                              parse_int_list(self.config.get('n_pcs_list', '0'))

        self.num_threads = int(self.config.get('num_threads', '4'))
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        # Where to save fitted models
        self.saved_models_dir = _resolve_outdir(
            base_outdir=self.model_output_dir,
            subdir="Saved_models",
            ensure_dir=True,
        )

    def _save_model(self, model_name: str, model_obj, file_index: str, n_pcs: int, n_snps: int, feature_list: List[str]):
        base = f"{model_name}_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs"
        feats_path = os.path.join(self.saved_models_dir, f"{base}.features.txt")
        with open(feats_path, "w") as f:
            for feat in feature_list:
                f.write(f"{feat}\n")

        if model_name == "MLP":
            path = os.path.join(self.saved_models_dir, f"{base}.pt")
            torch.save({"state_dict": model_obj.state_dict(), "feature_list": feature_list}, path)
            return path

        if model_name == "XGBoost":
            path = os.path.join(self.saved_models_dir, f"{base}.json")
            model_obj.save_model(path)
            return path

        # sklearn
        if joblib is not None:
            path = os.path.join(self.saved_models_dir, f"{base}.joblib")
            joblib.dump({"model": model_obj, "feature_list": feature_list}, path)
            return path

        import pickle
        path = os.path.join(self.saved_models_dir, f"{base}.pkl")
        with open(path, "wb") as f:
            pickle.dump({"model": model_obj, "feature_list": feature_list}, f)
        return path

    def run_machine_learning_and_save_results(self, args):
        """
        Adapted from the reference run_machine_learning_and_save_results.
        Adds: saving fitted model to disk.
        """
        model_name, train_X, train_y, test_X, test_y, file_index, n_pcs, n_snps, feature_list = args

        if patch_sklearn is not None:
            patch_sklearn(global_patch=True)

        # Initialize model
        if model_name == 'SVM':
            model = svm.SVC(kernel='rbf', probability=True)
        elif model_name == 'RandomForest':
            model = ensemble.RandomForestClassifier(n_estimators=100, max_depth=3)
        elif model_name == 'XGBoost':
            model = XGBClassifier(n_estimators=100, learning_rate=0.3)
        elif model_name == 'MLP':
            X_train_tensor = torch.tensor(train_X, dtype=torch.float32)
            X_test_tensor = torch.tensor(test_X, dtype=torch.float32)
            Y_train_tensor = torch.tensor(train_y, dtype=torch.float32).unsqueeze(1) if self.data_type == 'binary' else torch.tensor(train_y, dtype=torch.long)
            Y_test_tensor = torch.tensor(test_y, dtype=torch.float32).unsqueeze(1) if self.data_type == 'binary' else torch.tensor(test_y, dtype=torch.long)

            if self.data_type == 'binary':
                model = nnModel(num_i=len(feature_list), num_h=64, num_o=1).to(self.device)
                criterion = nn.BCEWithLogitsLoss()
            elif self.data_type == 'multiclass':
                model = nnModel(num_i=len(feature_list), num_h=64, num_o=2).to(self.device)
                criterion = nn.CrossEntropyLoss()
            else:
                raise ValueError("Unsupported data type for MLP.")
            optimizer = optim.Adam(model.parameters(), lr=0.01)
            batch_size = 32
            n_epochs = 100
            train_dataset = TensorDataset(X_train_tensor, Y_train_tensor)
            test_dataset = TensorDataset(X_test_tensor, Y_test_tensor)
            train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
            test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        else:
            raise ValueError(f"Model '{model_name}' is not supported.")

        loss_values = []
        if model_name == 'MLP':
            for epoch in range(n_epochs):
                model.train()
                correct = 0
                total = 0
                for X_batch, Y_batch in train_loader:
                    X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                    optimizer.zero_grad()
                    Y_pred = model(X_batch)
                    loss = criterion(Y_pred, Y_batch)
                    loss.backward()
                    optimizer.step()

                    if self.data_type == 'binary':
                        probabilities = torch.sigmoid(Y_pred)
                        predicted = (probabilities >= 0.5).float()
                    elif self.data_type == 'multiclass':
                        _, predicted = torch.max(Y_pred.data, 1)
                        Y_batch = Y_batch.squeeze()

                    total += Y_batch.size(0)
                    correct += (predicted.cpu() == Y_batch.cpu()).sum().item()

                loss_values.append(loss.item())
                accuracy_epoch = correct / total if total > 0 else 0.0
                print(f"[{model_name}] Epoch {epoch+1}/{n_epochs}, Loss: {loss.item():.4f}, Acc: {accuracy_epoch:.4f}")

            # Evaluate
            model.eval()
            with torch.no_grad():
                test_y_arr = []
                for X_batch, Y_batch in test_loader:
                    test_y_arr.extend(Y_batch.cpu().numpy())
                test_y_arr = np.array(test_y_arr).flatten()

                if self.data_type == "binary":
                    probs = torch.sigmoid(model(X_test_tensor.to(self.device))).detach().cpu().numpy().flatten()
                    probabilities = probs
                    predictions = (probs > 0.5).astype(int)
                else:
                    probs = nn.Softmax(dim=1)(model(X_test_tensor.to(self.device))).detach().cpu().numpy()
                    probabilities = probs
                    predictions = probs.argmax(axis=1)

                test_y = test_y_arr

        else:
            model.fit(train_X, train_y)
            if hasattr(model, "predict_proba"):
                probs = model.predict_proba(test_X)
                if self.data_type == 'binary':
                    probabilities = probs[:, 1]
                    predictions = (probabilities >= 0.5).astype(int)
                else:
                    probabilities = probs
                    predictions = np.argmax(probs, axis=1)
            else:
                probabilities = None
                predictions = model.predict(test_X)

        # Metrics
        accuracy = metrics.accuracy_score(test_y, predictions)
        mcc = metrics.matthews_corrcoef(test_y, predictions)
        f1 = metrics.f1_score(test_y, predictions)
        tn, fp, fn, tp = confusion_matrix(test_y, predictions).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        if unpatch_sklearn is not None:
            unpatch_sklearn()

        metrics_dict = {
            'model_name': model_name,
            'file_index': file_index,
            'n_pcs': n_pcs,
            'n_snps': n_snps,
            'accuracy': accuracy,
            'specificity': specificity,
            'mcc': mcc,
            'f1_score': f1,
        }
        metrics_df = pd.DataFrame([metrics_dict])

        results_file = os.path.join(self.model_output_dir, f'model_evaluation_results_{model_name}_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.csv')
        if not os.path.isfile(results_file):
            metrics_df.to_csv(results_file, index=False)
        else:
            metrics_df.to_csv(results_file, mode='a', header=False, index=False)

        np.save(os.path.join(self.model_output_dir, f'Probabilities_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs.npy'), probabilities)
        np.save(os.path.join(self.model_output_dir, f'Testy_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs.npy'), test_y)

        # NEW: save fitted model
        self._save_model(model_name, model, file_index, n_pcs, n_snps, feature_list)

        return metrics_dict

    def run_fold(self, file_fold_index: str, pcs: List[int], nsnps: List[int], models: List[str], max_workers: int):
        models_list = []
        for pc in pcs:
            for n_snp in nsnps:
                file_fold_npc_nsnp_index = f"{file_fold_index}_{pc}PCs_{n_snp}SNPs"
                current_dir = os.path.join(self.model_output_dir, file_fold_npc_nsnp_index)

                train_csv_path = os.path.join(current_dir, f'train_selected_snps_{file_fold_npc_nsnp_index}.csv')
                test_csv_path = os.path.join(current_dir, f'test_selected_snps_{file_fold_npc_nsnp_index}.csv')

                if not (os.path.exists(train_csv_path) and os.path.exists(test_csv_path)):
                    raise FileNotFoundError(f"Missing model input CSVs for {file_fold_npc_nsnp_index}")

                train_data = pd.read_csv(train_csv_path, index_col='taxa')
                test_data = pd.read_csv(test_csv_path, index_col='taxa')

                train_y = train_data['Phenotype'].values.astype(int)
                train_X = train_data.drop(columns=['Phenotype']).values
                test_y = test_data['Phenotype'].values.astype(int)
                test_X = test_data.drop(columns=['Phenotype']).values
                feature_list = train_data.drop(columns=['Phenotype']).columns.tolist()

                for model_name in models:
                    models_list.append((model_name, train_X, train_y, test_X, test_y, file_fold_index, pc, n_snp, feature_list))

        for item in models_list:
            self.run_machine_learning_and_save_results(item)

    def cat_result(self):
        result_files = [f for f in os.listdir(self.model_output_dir) if f.endswith('.csv')]
        all_data = []
        for file in result_files:
            file_path = os.path.join(self.model_output_dir, file)
            df = pd.read_csv(file_path)
            all_data.append(df)
        combined_df = pd.concat(all_data, ignore_index=True)
        sorted_df = combined_df.sort_values(by=['n_pcs', 'n_snps', 'model_name'])
        summarized_df = sorted_df.groupby(['n_pcs', 'n_snps', 'model_name'])[['accuracy', 'specificity', 'mcc', 'f1_score']].mean().reset_index()
        sorted_df.to_csv(os.path.join(self.output_dir, 'model_evaluation_results_all.csv'), index=False)
        summarized_df.to_csv(os.path.join(self.output_dir, 'model_evaluation_results_mean.csv'), index=False)
        print("Summary file saved")


class PredictionModelCVStep:
    def __init__(
        self,
        fold_configs: List[Path],
        max_workers: int,
        inner_model_workers: int,
        worker: Callable[[Dict[str, object]], str],
    ):
        self.fold_configs = fold_configs
        self.max_workers = max_workers
        self.inner_model_workers = inner_model_workers
        self.worker = worker

    def run(self) -> List[str]:
        finished: List[str] = []
        for cfg in self.fold_configs:
            finished.append(
                self.worker(
                    {
                        "fold_config": str(cfg),
                        "inner_model_workers": int(self.inner_model_workers),
                    }
                )
            )
        return finished
