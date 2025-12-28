#!/usr/bin/env python3
import subprocess,os,sys
import shap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearnex import patch_sklearn,unpatch_sklearn
from sklearn import metrics, ensemble, svm
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import confusion_matrix
from xgboost import XGBClassifier
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.init as init
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import roc_curve, auc
from torchmetrics import ROC
import argparse
from qmplot import manhattanplot, qqplot  # Ensure you have qmplot installed
from sklearn.preprocessing import StandardScaler
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple
from Utils.utils import read_config
import seaborn as sns
from imblearn.over_sampling import SMOTE



def parallel_process(items: List[tuple], process_func, max_workers: int = 8):
    """
    Perform parallel processing using ProcessPoolExecutor.
    Each item in `items` is a tuple of arguments for `process_func`.

    Args:
        items: List of tuples, where each tuple contains arguments for `process_func`.
        process_func: The function to execute in parallel.
        max_workers: Maximum number of parallel workers.
    """    
    print(f"Starting parallel processing with {max_workers} workers")
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_func, item) for item in items]
        for future in as_completed(futures):
            try:
                results.append(future.result())  # Collect results as they complete
            except Exception as e:
                print(f"Task failed with exception: {e}")

    print("Parallel processing complete")
    return results
    #with ProcessPoolExecutor(max_workers=max_workers) as executor:
    #    results = list(executor.map(process_func, items))
    #return results
class nnModel(nn.Module):
    def __init__(self, num_i, num_h, num_o):
        super(nnModel, self).__init__()
      
        self.fc1 = nn.Linear(num_i, num_h)
        self.fc2 = nn.Linear(num_h,num_h) #2 hidden layers
        self.fc3 = nn.Linear(num_h,num_o)

        init.xavier_normal_(self.fc1.weight)
        init.xavier_normal_(self.fc2.weight)
        init.xavier_normal_(self.fc3.weight)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)

        return x

class SNPPredictionModelTool:
    def __init__(self, config):
        """
        Initialize the model using the provided configuration dictionary.
        """
        # Use the provided configuration directly
        self.config = config
        
        # Set file paths
        self.species_name = self.config.get('species_name')
        self.phenotype_csv_path = self.config.get('phenotype_csv_path')
        self.vcf_file_path = self.config.get('vcf_file_path')
        self.chromosome_csv_path = self.config.get('chromosome_csv_path')
        self.phenotype_column = self.config.get('phenotype_column', 'Phenotype')
        # Get the working directory from the configuration file
        self.data_type = self.config.get('data_type', 'regression')
        self.oversampling = self.config.get('oversampling', 'True').lower() == 'true'

        base_output_dir = self.config.get('output_dir', 'results')
        prediction_output_dir = (self.config.get('prediction_output_dir') or '').strip()
        if prediction_output_dir:
            if os.path.isabs(prediction_output_dir):
                self.output_dir = prediction_output_dir
            else:
                self.output_dir = os.path.join(base_output_dir, prediction_output_dir)
        else:
            self.output_dir = base_output_dir
        self.blink_bin = self.config.get('blink_bin', 'blink1_5')
        blink_workers_raw = (self.config.get('blink_gwas_workers') or '').strip()
        try:
            self.blink_gwas_workers = int(blink_workers_raw) if blink_workers_raw else 1
        except ValueError:
            self.blink_gwas_workers = 1
        if self.blink_gwas_workers < 1:
            self.blink_gwas_workers = 1
        # Ensure the working directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load chromosome recode mapping
        self.load_chromosome_recode()
        
        # PCA and GWAS settings
        self.use_pcs = self.config.get('use_pcs', 'True').lower() == 'true'
        self.n_pcs_list = [int(x) for x in self.config.get('n_pcs_list', '0').split(',')]
        self.top_n_snps_list = [int(x) for x in self.config.get('top_n_snps_list', '10').split(',')]
        self.models = self.config.get('models', 'SVM,RandomForest,XGBoost').split(',')
        self.num_shuffles = int(self.config.get('num_shuffles', '10'))
        threads_raw = (self.config.get('num_threads') or self.config.get('threads') or '4').strip()
        try:
            self.num_threads = int(threads_raw)
        except ValueError:
            self.num_threads = 4
        self.num_splits = int(self.config.get('num_splits', '5'))
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        

        # Create a unique temporary directory within the working directory
        self._create_directories()

    def _create_directories(self):
        """
        Create necessary directories and subdirectories.
        """
        # Create the main directories
        os.makedirs(self.output_dir, exist_ok=True)

        # Create subdirectories
        os.makedirs(os.path.join(self.output_dir, 'PCA_Scree_Plot'), exist_ok=True)
        self.pca_scree_output_dir = os.path.join(self.output_dir, 'PCA_Scree_Plot')

        os.makedirs(os.path.join(self.output_dir, 'GWAS_dir'), exist_ok=True)
        self.gwas_output_dir = os.path.join(self.output_dir, 'GWAS_dir')

        os.makedirs(os.path.join(self.output_dir, 'Shap_dir'), exist_ok=True)
        self.shap_output_dir = os.path.join(self.output_dir, 'Shap_dir')

        os.makedirs(os.path.join(self.output_dir, 'AUC_dir'), exist_ok=True)
        self.auc_plot_output_dir = os.path.join(self.output_dir, 'AUC_dir')

        os.makedirs(os.path.join(self.output_dir, 'Model_result'), exist_ok=True)
        self.model_output_dir = os.path.join(self.output_dir, 'Model_result')

        print(f"Directories created under '{self.output_dir}'.")

    
    def load_chromosome_recode(self):
        """
        Load the chromosome recode CSV file into mappings.
        """
        self.chromosome_mapping = pd.read_csv(self.chromosome_csv_path)
        # Create mapping dictionaries
        self.chr_name_to_recoded = dict(zip(self.chromosome_mapping['OriginalName'], self.chromosome_mapping['RecodedNumber']))
        self.recoded_to_label = dict(zip(self.chromosome_mapping['RecodedNumber'].astype(str), self.chromosome_mapping['Label']))

    def vcf_to_dataframe(self, vcf_file):
        """
        Convert a VCF file to a pandas DataFrame.
        """
        # Extract header
        header = self.get_vcf_header(vcf_file)
        # Read VCF using pandas
        vcf_reader = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, dtype={0: str})
        vcf_reader.columns = header
        # Identify sample columns (columns from index 9 onwards)
        sample_columns = header[9:]
        # Extract genotype data for samples
        genotype_df = vcf_reader[['ID'] + sample_columns]
        # Transpose the DataFrame so that samples are rows and SNPs are columns
        genotype_df = genotype_df.transpose()
        # Set sample IDs as index
        genotype_df.columns = genotype_df.iloc[0]
        genotype_df = genotype_df[1:]
        return genotype_df
    
    def get_vcf_header(self, vcf_file):
        """
        Extract the header from a VCF file.
        """
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    return line.strip().split('\t')
        raise ValueError("VCF header not found.")
    
    def convert_genotypes_to_numeric(self, genotype_df):
        """
        Convert genotype strings to numeric values and also fill in NA.

        Parameters:
        - genotype_df: pandas DataFrame with genotype strings.

        Returns:
        - genotype_numeric_df: pandas DataFrame with numeric genotypes.
        """

        genotype_numeric_df = genotype_df.apply(lambda col: col.map(lambda x: 1 if x == '0/0' or x == '0|0'
                                                               else 2 if x == '0/1' or x == '0|1'
                                                               else 3 if x == '1/1' or x == '1|1'
                                                               else np.nan if x == './.' 
                                                               else x))
        # Impute missing values
        # Convert all values to float
        genotype_numeric_df = genotype_numeric_df.astype(float)
        imputer = SimpleImputer(strategy='median')
        genotype_imputed = pd.DataFrame(imputer.fit_transform(genotype_numeric_df), 
                                        columns=genotype_numeric_df.columns, 
                                        index=genotype_numeric_df.index)
        return genotype_imputed
    
    def perform_pca(self, genotype_df, fold_index, n_components=15):
        """
        Perform PCA on the genotype data.
        """
        # Standardize the data
        scaler = StandardScaler()
        genotype_scaled = scaler.fit_transform(genotype_df)
        # Perform PCA
        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(genotype_scaled)
        explained_variance_ratio = pca.explained_variance_ratio_
        # Create a DataFrame with principal components
        pc_columns = [f'PC{i+1}' for i in range(n_components)]
        principal_df = pd.DataFrame(data=principal_components, columns=pc_columns)
        principal_df.insert(0, 'taxa', genotype_df.index)  # (index, column_name, values)
        
        return principal_df, explained_variance_ratio
    
    def select_top_snps(self, gwas_results_file, top_n_snps=10):
        """
        Select top N SNPs based on p-value.
        """
        gwas_df = pd.read_csv(gwas_results_file, sep='\t')
        gwas_df = gwas_df[gwas_df['p_value'] > 0]
        gwas_df['-logp'] = -np.log10(gwas_df['p_value'])
        sorted_df = gwas_df.sort_values(by='-logp', ascending=False)
        top_n_unique_values = sorted_df['-logp'].unique()[:top_n_snps]
        top_df = sorted_df[sorted_df['-logp'].isin(top_n_unique_values)]
        top_snps = top_df['taxa'].tolist()
        return top_snps
    
    def blink_file_generator(self, train_df, file_fold_npc_index, principal_df, n_pcs):
        """
        Generate input files for BLINK GWAS analysis.
        """
        train_samples = train_df.index.tolist()
        train_samples_str = ','.join(train_samples)

        # File directory
        gwas_output_dir = os.path.join(self.output_dir, 'GWAS_dir')
        os.makedirs(os.path.join(gwas_output_dir, f'{file_fold_npc_index}'), exist_ok=True)
        current_dir = os.path.join(gwas_output_dir, f'{file_fold_npc_index}')

        # Generate phenotype file
        phenotype_df = pd.read_csv(self.phenotype_csv_path,sep=None, engine='python')
        phenotype_df = phenotype_df[['taxa', self.phenotype_column]]
        filtered_phenotype_df = phenotype_df[phenotype_df['taxa'].isin(train_samples)]
        phenotype_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.txt')
        filtered_phenotype_df.to_csv(phenotype_file, sep='\t', index=False)
    
        # Generate genotype file
        vcf_file = self.vcf_file_path
        output_vcf = os.path.join(current_dir, f'train_{file_fold_npc_index}.vcf')
        mapping_file = os.path.join('chr_mapping.txt')
        # Create chromosome mapping file
        # Create chromosome mapping file
        if not os.path.exists(mapping_file):
            self.chromosome_mapping[['OriginalName', 'RecodedNumber']].to_csv(
                                    mapping_file, sep='\t', header=False, index=False
        )

        # Get the list of original chromosome names to keep
        original_chromosomes_to_keep = self.chromosome_mapping['OriginalName'].tolist()
        regions_str = ','.join(map(str, original_chromosomes_to_keep))

        # Create a temporary compressed file
        compressed_vcf = vcf_file + ".gz"
        if not os.path.exists(compressed_vcf):
            subprocess.run(f"bgzip -c {vcf_file} > {compressed_vcf}", shell=True, check=True)
            # Index the compressed file
            subprocess.run(f"bcftools index {compressed_vcf}", shell=True, check=True)
    
    
        # Build the bcftools command
        command = (
            f'bcftools view -s "{train_samples_str}" -r "{regions_str}" "{compressed_vcf}" | \
              bcftools annotate --rename-chrs "{mapping_file}" -Oz -o "{output_vcf}"'
        )

        # Run the command
        subprocess.run(command, shell=True, check=True, executable="/bin/bash")

        # Generate covariate file (if PCs are used)
        if self.use_pcs and n_pcs > 0:
            principal_df = principal_df.iloc[:, :(n_pcs+1)]
            filtered_principal_df = principal_df[principal_df['taxa'].isin(train_samples)]
            principal_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.cov')
            filtered_principal_df.to_csv(principal_file, sep='\t', index=False)
            # Save explained variance ratios for scree plot

    def run_gwas(self, file_fold_npc_index):
        """
        Run GWAS analysis using BLINK.
        """
        # File directory
        current_dir = os.path.join(self.gwas_output_dir, f'{file_fold_npc_index}')
        os.makedirs(current_dir, exist_ok=True)
        tag = f"train_{file_fold_npc_index}"
        # Remove '> blink_run_{tag}.log' from the command
        blink_command = (
            f'{self.blink_bin} '
            f'--gwas '
            f'--file {tag} '
            f'--vcf '
            f'--out {tag}'
        )
        # Open the log file in write mode and redirect both stdout and stderr
        log_path = os.path.join(current_dir, f'blink_run_{tag}.log')
        with open(log_path, 'w') as log_file:
            subprocess.run(
                blink_command,
                shell=True,
                check=True,
                cwd=current_dir,
                stdout=log_file,            # Redirect STDOUT to the log file
                stderr=subprocess.STDOUT    # Redirect STDERR to the same log file
            )


    def select_markers_and_prepare_data_for_model_construction(self, file_fold_npc_index, file_fold_npc_nsnp_index, test_df, top_n_snps=10):
        """
        Select markers and prepare train and test datasets, incorporating phenotype data.
        """
        current_dir = os.path.join(self.gwas_output_dir, file_fold_npc_index)
        p_value_file = os.path.join(current_dir, f'train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt')
        vcf_file = os.path.join(current_dir, f'train_{file_fold_npc_index}.vcf')

        # Select top SNPs
        top_snps = self.select_top_snps(p_value_file, top_n_snps)

        # Read phenotype data
        phenotype_df = pd.read_csv(self.phenotype_csv_path,sep=None, engine='python')
        phenotype_df = phenotype_df[['taxa', self.phenotype_column]]
        phenotype_df.set_index('taxa', inplace=True)

        # Prepare train data
        train_genotype_df = self.vcf_to_dataframe(vcf_file)
        train_selected_snps = train_genotype_df[top_snps]
        train_selected_snps = self.convert_genotypes_to_numeric(train_selected_snps)
        
        # Prepare test data
        test_selected_snps = test_df[top_snps]
        test_selected_snps = self.convert_genotypes_to_numeric(test_selected_snps)

        # Merge phenotype data
        train_selected_snps = train_selected_snps.join(phenotype_df, how='inner')
        train_selected_snps.rename(columns={self.phenotype_column: 'Phenotype'}, inplace=True)

        test_selected_snps = test_selected_snps.join(phenotype_df, how='inner')
        test_selected_snps.rename(columns={self.phenotype_column: 'Phenotype'}, inplace=True)

        #Save directory
        save_dir = os.path.join(self.model_output_dir, f'{file_fold_npc_nsnp_index}')
        os.makedirs(save_dir, exist_ok=True)
        # Save train data
        train_csv_path = os.path.join(save_dir, f'train_selected_snps_{file_fold_npc_nsnp_index}.csv')
        train_selected_snps.to_csv(train_csv_path, index_label='taxa')

        # Save test data
        test_csv_path = os.path.join(save_dir, f'test_selected_snps_{file_fold_npc_nsnp_index}.csv')
        test_selected_snps.to_csv(test_csv_path, index_label='taxa')

        return top_snps
    

    def run_machine_learning_and_save_results(self, args):
        """
        Run machine learning models, save SHAP plots, and save results to a CSV file.

        Args:
            args: Tuple containing:
                model_name (str): Name of the machine learning model.
                train_X, train_y: Training data and labels.
                test_X, test_y: Test data and labels.
                file_index (str): File identifier.
                n_pcs (int): Number of principal components.
                n_snps (int): Number of SNPs.
                feature_list (list): List of features.
        """
        model_name, train_X, train_y, test_X, test_y, file_index, n_pcs, n_snps, feature_list = args
        # print(f"Task {args[0]} (Model: {args[0]}) is running in Process ID: {os.getpid()};file_index:{file_index};n_pcs:{n_pcs};n_snps:{n_snps}")
        # Initialize model
        patch_sklearn(global_patch=True)
        if model_name == 'SVM':
            model = svm.SVC(kernel='rbf', probability=True)  
        elif model_name == 'RandomForest':
            model = ensemble.RandomForestClassifier(n_estimators=100, max_depth=3)
        elif model_name == 'XGBoost':
            model = XGBClassifier(n_estimators=100, learning_rate=0.3)
        elif model_name == 'MLP':
            X_train_tensor = torch.tensor(train_X, dtype = torch.float32)
            X_test_tensor = torch.tensor(test_X, dtype = torch.float32)
            Y_train_tensor = torch.tensor(train_y, dtype=torch.float32).unsqueeze(1) if self.data_type == 'binary' else torch.tensor(train_y, dtype=torch.long)
            Y_test_tensor = torch.tensor(test_y, dtype=torch.float32).unsqueeze(1) if self.data_type == 'binary' else torch.tensor(test_y, dtype=torch.long)
            
            if self.data_type == 'binary':
                model = nnModel(num_i=len(feature_list), num_h=64, num_o=1).to(self.device)  # Single output
                criterion = nn.BCEWithLogitsLoss()
            elif self.data_type == 'multiclass':
                model = nnModel(num_i=len(feature_list), num_h=64, num_o=2).to(self.device)  # Multi-class output
                criterion = nn.CrossEntropyLoss()
            else:
                raise ValueError("Unsupported data type for MLP.")
            optimizer = optim.Adam(model.parameters(), lr = 0.01)
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
                        _,predicted = torch.max(Y_pred.data, 1)
                        Y_batch = Y_batch.squeeze()

                    total += Y_batch.size(0)
                    correct += (predicted.cpu() == Y_batch.cpu()).sum().item()

                    
                loss_values.append(loss.item())
                accuracy = correct/total
                print(f"Epoch {epoch+1}/{n_epochs}, Loss: {loss.item():.4f}, Acc: {accuracy:.4f}")

            # Plot loss over epochs
            plt.figure()
            plt.plot(range(1, n_epochs + 1), loss_values, marker='o')
            plt.title(f'{model_name} Training Loss over Epochs (Fold {file_index}, PC{n_pcs}, {n_snps} SNPs)')
            plt.xlabel('Epoch')
            plt.ylabel('Loss')
            loss_plot_path = os.path.join(self.model_output_dir, f'{model_name}_Loss_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.png')
            plt.savefig(loss_plot_path)
            plt.close()

            model.eval()

            with torch.no_grad():
                predictions = []
                test_y = []

                for X_batch, Y_batch in test_loader:
                    X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                    Y_pred = model(X_batch)

                    _,predicted = torch.max(Y_pred.data, 1)

                    predictions.extend(predicted.numpy())
                    test_y.extend(Y_batch.numpy())

                if self.data_type == "binary":
                    Y_batch = Y_batch.unsqueeze(1) 
                    probabilities = torch.sigmoid(model(X_test_tensor)).detach().numpy().flatten()
                    predictions = (probabilities > 0.5).astype(int)
                    
                    
                elif self.data_type == "multiclass":
                    probabilities = nn.Softmax(dim=1)(model(X_test_tensor)).detach().numpy()
                    predictions = probabilities.argmax(axis=1)
                
                test_y = np.array(test_y).flatten()

        else:
            # Fit model
            model.fit(train_X, train_y)

            if hasattr(model, "predict_proba"):
                probabilities = model.predict_proba(test_X)  # Shape (n_samples, n_classes)
                # predict_proba=True => largely increase the training time ~4-5 times
                if self.data_type == 'binary':
                    probabilities_class1 = probabilities[:, 1]
                    predictions = (probabilities_class1 >= 0.5).astype(int)
                    probabilities = probabilities_class1
                    
                else:
                    predictions = np.argmax(probabilities, axis=1)
            else:
                probabilities = None
                predictions = model.predict(test_X)

        
        # Metrics
        accuracy = metrics.accuracy_score(test_y, predictions)
        mcc = metrics.matthews_corrcoef(test_y, predictions)
        f1 = metrics.f1_score(test_y, predictions)
        tn, fp, fn, tp = confusion_matrix(test_y, predictions).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        # SHAP values
        if model_name in ['RandomForest', 'XGBoost']:
            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(test_X)
            # For binary classification, shap_values can be an array or list
            if model_name == 'RandomForest':
                shap_values_to_plot = shap_values[:,:,1]  # Assuming class 1 is the positive class
            else:
                shap_values_to_plot = shap_values

        elif model_name == "MLP":
            # Ensure tensors are on the CPU for SHAP
            
            model_cpu = model.cpu()
            X_train_cpu = X_train_tensor.cpu()
            X_test_cpu = X_test_tensor.cpu()
            
            explainer = shap.GradientExplainer(model_cpu, X_train_cpu)
            shap_values = explainer.shap_values(X_test_cpu)
            shap_values_to_plot = np.squeeze(shap_values)

            # Move model back to original device
            model.to(self.device)

        else:
            # For models without native SHAP support, use KernelExplainer
            explainer = shap.KernelExplainer(model.predict_proba, train_X)
            shap_values = explainer.shap_values(test_X)
            shap_values_to_plot = shap_values[:,:,1]  # For binary classification

        unpatch_sklearn()

        # Save the DataFrame to a CSV file
        df_shap_values = pd.DataFrame(shap_values_to_plot, columns=feature_list) 
        df_feature_importance = pd.DataFrame(columns=['feature','importance'])
        for col in df_shap_values.columns:
            importance = df_shap_values[col].abs().mean()
            df_feature_importance.loc[len(df_feature_importance)] = [col,importance]
            df_feature_importance = df_feature_importance.sort_values('importance',ascending=False)   
        
        df_feature_importance_file = os.path.join(self.model_output_dir, f'ShapTop10_{model_name}_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.csv')
        df_feature_importance.to_csv(df_feature_importance_file, index=False)

        # Plot SHAP summary plot
        plt.figure()
        shap.summary_plot(shap_values_to_plot, test_X, feature_names=feature_list, show=False)
        plt.title(f'{model_name} SHAP Summary Plot (Fold {file_index},PC{n_pcs},{n_snps}SNPs)')
        shap_plot = os.path.join(self.shap_output_dir, f'{model_name}_SHAP_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.png')
        plt.savefig(shap_plot)
        plt.close()
        
        
        # Prepare metrics dictionary
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

        # Create DataFrame
        metrics_df = pd.DataFrame([metrics_dict])

        # Save metrics to CSV file
        results_file = os.path.join(self.model_output_dir, f'model_evaluation_results_{model_name}_Fold_{file_index}_PC{n_pcs}_{n_snps}SNPs.csv')
        if not os.path.isfile(results_file):
            # If file does not exist, write header
            metrics_df.to_csv(results_file, index=False)
        else:
            # Append to existing file
            metrics_df.to_csv(results_file, mode='a', header=False, index=False)

        # Optionally, save probabilities and test_y to separate files for future use
        #print(f'Probabilities_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs', probabilities.shape)
        #print(f'Testy_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs', test_y.shape)
        np.save(os.path.join(self.model_output_dir, f'Probabilities_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs.npy'), probabilities)
        np.save(os.path.join(self.model_output_dir, f'Testy_{model_name}_Fold{file_index}_PC{n_pcs}_{n_snps}SNPs.npy'), test_y)

    def cat_result(self):
        """
        Concatenate and summarize all result CSV files in the specified directory.
        Compute the mean of specific columns for rows matching the given condition.

        Parameters:
        - output_dir (str): Directory containing the result CSV files.
        - condition_column (str): Column name to apply the condition on.
        - condition_value (any): Value of the condition to filter rows.
        - mean_columns (list of str): Columns to compute the mean for.
        - summary_file (str): Path to save the clean summarized output file.

        Returns:
        - pd.DataFrame: Summarized DataFrame.
        """
        # Find all result files in the directory
        result_files = [f for f in os.listdir(self.model_output_dir) if f.endswith('.csv')]
        all_data = []

        # Read and concatenate all files
        for file in result_files:
            file_path = os.path.join(self.model_output_dir, file)
            df = pd.read_csv(file_path)
            all_data.append(df)

        # Combine all dataframes into a single dataframe
        combined_df = pd.concat(all_data, ignore_index=True)

        # Sort the data by the specified column
        sorted_df = combined_df.sort_values(by=['n_pcs', 'n_snps', 'model_name'])

        # Group by the specified column and compute mean for the desired columns
        summarized_df = sorted_df.groupby(['n_pcs', 'n_snps', 'model_name'])[['accuracy', 'specificity', 'mcc', 'f1_score']].mean().reset_index()

        print(summarized_df)
        # Save the summarized data to a file
        sorted_df.to_csv(os.path.join(self.output_dir, 'model_evaluation_results_all.csv'), index=False)
        summarized_df.to_csv(os.path.join(self.output_dir, 'model_evaluation_results_mean.csv'), index=False)

        print(f"Summary file saved")


    def plot_performance_comparison(self):
        """
        Generate comparison plots for n_snps vs. models and n_pcs vs. models in each fold,
        reading data from saved .npy files for probabilities and true labels.

        The plots show the performance (AUC) of different models for varying numbers of SNPs and PCs.

        Outputs:
        - AUC comparison plots for n_snps vs. models and n_pcs vs. models.
        - A summary CSV file with the best AUC values per model and configuration.
        """
        # Find all saved .npy files in the output directory
        probabilities_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Probabilities")]
        testy_files = [f for f in os.listdir(self.model_output_dir) if f.startswith("Testy")]

        # Ensure probabilities and test_y files match
        probabilities_files.sort()
        testy_files.sort()

        auc_scores = []  # To store AUC scores for all configurations

        # Loop through each saved result
        for prob_file, testy_file in zip(probabilities_files, testy_files):
            # Extract metadata from the file names
            meta = prob_file.split('_')
            model_name = meta[1]
            fold = next((frag.replace('Fold', '') for frag in meta if frag.startswith('Fold')), None)
            n_pcs = int(next((frag.replace('PC', '') for frag in meta if 'PC' in frag), None))
            n_snps = int(next((frag.replace('SNPs.npy', '') for frag in meta if 'SNPs.npy' in frag), None))

            # Load the probabilities and test_y arrays
            y_scores = np.load(os.path.join(self.model_output_dir, prob_file))
            test_y = np.load(os.path.join(self.model_output_dir, testy_file))

            # Calculate ROC curve and AUC
            fpr, tpr, _ = roc_curve(test_y, y_scores)

            # Store AUC score with metadata
            auc_scores.append({
                'model_name': model_name,
                'fold': fold,
                'n_pcs': n_pcs,
                'n_snps': n_snps,
                'fpr': fpr,
                'tpr': tpr
            })

        # Convert AUC scores to DataFrame for easier analysis
        auc_df = pd.DataFrame(auc_scores)

        #roc_auc = auc(auc_df['fpr'], auc_df['tpr'])
        auc_df['roc_auc'] = auc_df.apply(
            lambda row: auc(row['fpr'], row['tpr']) if isinstance(row['fpr'], (list, np.ndarray)) and 
            isinstance(row['tpr'], (list, np.ndarray)) and len(row['fpr']) == len(row['tpr']) else np.nan,
            axis=1
            )
        #auc_df['roc_auc'] = roc_auc
        
        # Define consistent color and line style mappings
        model_colors = sns.color_palette("tab10", n_colors=len(auc_df['model_name'].unique()))
        model_color_map = {model: model_colors[i] for i, model in enumerate(auc_df['model_name'].unique())}

        line_styles = ['solid', 'dashed', 'dotted', 'dashdot',(0, (1, 1)), # Densely dotted
                        (0, (5, 5)),            # Loosely dashed
                        (0, (3, 5, 1, 5))]      # Dash-dot-dash]

        # Plot n_snps vs. models
        for fold in auc_df['fold'].unique():
            for n_snp in auc_df['n_snps'].unique():
                plt.figure(figsize=(6, 6))
                fold_nsnp_data = auc_df[(auc_df['fold'] == fold) & (auc_df['n_snps'] == n_snp) ]
                sorted_data = fold_nsnp_data.sort_values(by=['model_name', 'n_pcs'])
                plt.plot([0, 1], [0, 1], 'k--', label='Random (y=x)')
                for _, row in sorted_data.iterrows():  # Iterate over rows
                    model = row['model_name']
                    fpr = row['fpr']  # Extract list
                    tpr = row['tpr']  # Extract list
                    npc_or_nsnps_styles = {n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df['n_pcs'].unique()))}

                    # Use consistent color for model and line style for n_snps
                    plt.plot(fpr,tpr,label=f"{model}, {row['n_pcs']}PCs (AUC = {row['roc_auc']:.2f})",
                            color=model_color_map[model],
                            linestyle=npc_or_nsnps_styles[row['n_pcs']]
                            )

                plt.title(f'ROC curve (Fold {fold}, {n_snp}SNPs): PCs vs. Models of {self.phenotype_column} in {self.species_name}')
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.grid()

                # Place the legend outside the plot
                plt.legend(
                    loc='upper left',
                    bbox_to_anchor=(1, 1),
                    title='Legend',
                    fontsize='small',
                    markerscale=0,  # Scale marker size in legend
                    handletextpad=0.8 # Padding between marker and text
                )
                output_path = os.path.join(self.auc_plot_output_dir, f'auc_comparison_fold{fold}_{n_snp}SNPs.png')
                plt.savefig(output_path, bbox_inches='tight')
                plt.close()
                print(f'Plot saved: {output_path}')

            for n_pc in auc_df['n_pcs'].unique():
                plt.figure(figsize=(6, 6))
                fold_npc_data = auc_df[(auc_df['fold'] == fold) & (auc_df['n_pcs'] == n_pc) ]

                # Add diagonal line
                plt.plot([0, 1], [0, 1], 'k--', label='Random (y=x)')

                for _, row in fold_npc_data.iterrows():  # Iterate over rows
                    model = row['model_name']
                    fpr = row['fpr']  # Extract list
                    tpr = row['tpr']  # Extract list

                    npc_or_nsnps_styles = {n: line_styles[i % len(line_styles)] for i, n in enumerate(sorted(auc_df['n_snps'].unique()))}
                    # Use consistent color for model and line style for n_snps
                    plt.plot(fpr,tpr,label=f"{model}, {row['n_snps']}SNPs (AUC = {row['roc_auc']:.2f})",
                            color=model_color_map[model],
                            linestyle=npc_or_nsnps_styles[row['n_snps']]
                            )

                plt.title(f'ROC curve (Fold {fold}, {n_pc}PCs): SNPs vs. Models of {self.phenotype_column} in {self.species_name}')
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.grid()

                # Place the legend outside the plot
                plt.legend(
                    loc='upper left',
                    bbox_to_anchor=(1, 1),
                    title='Legend',
                    fontsize='small',
                    markerscale=0,  # Scale marker size in legend
                    handletextpad=0.8 # Padding between marker and text
                )
                output_path = os.path.join(self.auc_plot_output_dir, f'auc_comparison_fold{fold}_{n_pc}PCs.png')
                plt.savefig(output_path, bbox_inches='tight')
                plt.close()
                print(f'Plot saved: {output_path}')

        # Save AUC summary to CSV
        auc_df[['fpr', 'tpr']] = auc_df[['fpr', 'tpr']].applymap(lambda x: f'{x}')
        summary_csv = os.path.join(self.auc_plot_output_dir, 'auc_summary.csv')
        auc_df.to_csv(summary_csv, index=False)
        print(f'Summary CSV saved: {summary_csv}')

    def pca_plotter(self, pca_df_by_fold):
        """
        Plot the first three principal components from the PCA results across folds.
        """
        # Set markers for the groups
        markers = ['o', 's', '^', 'D', 'v', '*', '+', 'x']  # Extend markers if needed

        # Create a figure with three subplots for PC1 vs PC2, PC2 vs PC3, PC1 vs PC3
        fig, axes = plt.subplots(1, 3, figsize=(24, 8))

        # Titles for each plot
        titles = ['PC1 vs PC2', 'PC2 vs PC3', 'PC1 vs PC3']

        # Define colors for each group
        unique_groups = None  # Will be determined based on phenotype data

        for fold_index, df in pca_df_by_fold.items():
            # Read phenotype data to get group labels
            phenotype_df = pd.read_csv(self.phenotype_csv_path, sep=None, engine='python')
            phenotype_df = phenotype_df[['taxa', self.phenotype_column]]
            df = df.merge(phenotype_df, on='taxa')

            # Determine unique groups if not already determined
            if unique_groups is None:
                unique_groups = sorted(df[self.phenotype_column].unique())

            # Plot PC1 vs PC2, PC2 vs PC3, and PC1 vs PC3 for each group
            for group_value, marker in zip(unique_groups, markers):
                group_df = df[df[self.phenotype_column] == group_value]

                # PC1 vs PC2
                axes[0].scatter(group_df['PC1'], group_df['PC2'], s=50, 
                                label=f'Fold {fold_index} - Group {group_value}', marker=marker)

                # PC2 vs PC3
                axes[1].scatter(group_df['PC2'], group_df['PC3'], s=50, 
                                label=f'Fold {fold_index} - Group {group_value}', marker=marker)

                # PC1 vs PC3
                axes[2].scatter(group_df['PC1'], group_df['PC3'], s=50, 
                                label=f'Fold {fold_index} - Group {group_value}', marker=marker)

        # Set titles and labels for each subplot
        axes[0].set_title(titles[0], fontsize=16)
        axes[0].set_xlabel('Principal Component 1', fontsize=14)
        axes[0].set_ylabel('Principal Component 2', fontsize=14)

        axes[1].set_title(titles[1], fontsize=16)
        axes[1].set_xlabel('Principal Component 2', fontsize=14)
        axes[1].set_ylabel('Principal Component 3', fontsize=14)

        axes[2].set_title(titles[2], fontsize=16)
        axes[2].set_xlabel('Principal Component 1', fontsize=14)
        axes[2].set_ylabel('Principal Component 3', fontsize=14)

        # Add legends to each subplot and position them outside the plot
        for ax in axes:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Legend')  # Move legend outside
            ax.grid(True)

        # Save the plot
        plt.tight_layout()
        plot_path = os.path.join(self.pca_scree_output_dir, f'pca_plots_folds_{self.phenotype_column}_{self.species_name}.png')
        plt.savefig(plot_path)
        plt.close()
    
    def combined_scree_plot(self, explained_variances_by_fold):
        """
        Plot the scree plot of explained variance for the first 15 principal components across folds.
    
        Parameters:
        - explained_variances_by_fold (dict): Dictionary where keys are fold indices and values are explained variance ratios.
        """
        # Create a figure for the scree plot
        fig, ax = plt.subplots(figsize=(8, 6))

        # Define colors for each fold
        num_folds = len(explained_variances_by_fold.keys())
        colors = plt.cm.rainbow(np.linspace(0, 1, num_folds))

        # Map fold indices to sequential indices
        fold_index_map = {fold: idx for idx, fold in enumerate(explained_variances_by_fold.keys())}

        # Loop through the dictionary
        for fold_index, explained_variance_ratio in explained_variances_by_fold.items():
            # Map fold_index to sequential color index
            color_index = fold_index_map[fold_index]

            # Extract explained variance ratios for the first 15 PCs
            cumulative_variance = np.cumsum(explained_variance_ratio[:15]) * 100

            # Plot explained variance for PC1-PC15 for each fold in the same plot
            ax.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, marker='o',
                color=colors[color_index], label=f'Fold {fold_index}')

        # Adding titles and labels
        ax.set_xticks(range(1, 16))
        ax.set_title(f'Scree Plot of Explained Variance for PCA - {self.phenotype_column} in {self.species_name}', fontsize=16)
        ax.set_xlabel('Principal Components', fontsize=14)
        ax.set_ylabel('Cumulative Explained Variance Ratio (%)', fontsize=14)

        # Display the legend for the fold colors
        ax.legend(loc='lower right', title='Folds')
        plt.grid(True)
        plt.tight_layout()

        # Save the plot to a file
        plot_path = os.path.join(self.pca_scree_output_dir, f'scree_plot_folds_{self.phenotype_column}_{self.species_name}.png')
        plt.savefig(plot_path)
        plt.close()
    
    def plot_manhattan_and_qq(self, file_fold_npc_index):
        """
        Generate Manhattan and QQ plots for GWAS results.
        """

        current_dir = os.path.join(self.gwas_output_dir, file_fold_npc_index)
        os.makedirs(current_dir, exist_ok=True)
        
        p_value_file = os.path.join(current_dir, f'train_{file_fold_npc_index}_{self.phenotype_column}_GWAS_result.txt')
        
        df = pd.read_table(p_value_file, sep="\t")

        # Sort DataFrame
        df = df.sort_values(by=['chr', 'pos'])

        # Rename columns
        df = df.rename(columns={'chr':'#CHROM','pos': 'POS','p_value': 'P','taxa': 'ID' })

        # Map chromosome names using recode_dict
        recode_dict = dict(zip(self.chromosome_mapping['OriginalName'], self.chromosome_mapping['Label']))
        df['#CHROM'] = df['#CHROM'].map(recode_dict)

        # Remove rows with P-value of zero
        df = df[df["P"] != 0]

        # Generate the Manhattan plot
        fig, ax = plt.subplots(figsize=(12, 6))
        ax = manhattanplot(
            data=df,
            sign_line_cols=["#D62728", "#2CA02C"],
            hline_kws={"linestyle": "--", "lw": 1.3},
            xticklabel_kws={"rotation": 45},
            sign_marker_p=0.05/df.shape[0],
            sign_marker_color="r",
            title=f"Manhattan Plot - Fold {file_fold_npc_index} of {self.phenotype_column} in {self.species_name}",
            xlabel="Chromosome",
            ylabel=r"$-log_{10}{(P)}$",
            ax=ax
        )
        output_manplot_path = os.path.join(current_dir, f'manhattan_plot_{file_fold_npc_index}.png')
        plt.savefig(output_manplot_path)
        plt.close()

        # Generate the QQ plot
        plt.figure(figsize=(6, 6))
        qqplot(data=df["P"])
        plt.title(f"QQ Plot - Fold {file_fold_npc_index} of {self.phenotype_column} in {self.species_name}")
        output_qqplot_path = os.path.join(current_dir, f'qq_plot_shuffle-fold_{file_fold_npc_index}.png')
        plt.savefig(output_qqplot_path)
        plt.close()

def run_snp_prediction_model():
    """
    Run SNP Prediction Model with x-fold cross-validation, n_snps, n_pcs, models.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="SNP Prediction Model")
    parser.add_argument("config_file", help="Path to the configuration file")
    parser.add_argument(
        "--stage",
        choices=["all", "gwas", "model"],
        default="all",
        help="Run BLINK GWAS/PCA, model construction, or both.",
    )
    args = parser.parse_args()
    print("Path to config file:", args.config_file)

    # Read configuration file
    config = read_config(args.config_file)

    # Initialize the model with the configuration
    Tools = SNPPredictionModelTool(config)

    stage = args.stage.lower()
    run_gwas = stage in {"all", "gwas"}
    run_model = stage in {"all", "model"}

    # Extract parameters
    eva_file_path = os.path.join(Tools.output_dir, 'model_evaluation_results_mean.csv')

    if os.path.isfile(eva_file_path) and stage in {"all", "model"}:
        print("Evaluation results already exist. Run the ROC evaluation script to plot performance.")
        print("Finsihed")
        return

    models = Tools.models if run_model else []

    # Load data
    genotype_df = Tools.vcf_to_dataframe(Tools.vcf_file_path)
    phenotype_df = pd.read_csv(Tools.phenotype_csv_path)
    phenotype_df.set_index('taxa', inplace=True)
    y = phenotype_df.loc[genotype_df.index, Tools.phenotype_column].values

    # Initialize StratifiedKFold for standard x-fold cross-validation
    skf = StratifiedKFold(n_splits=Tools.num_splits, shuffle=True, random_state=2024)
    splits = list(skf.split(genotype_df, y))
    explained_variances_by_fold = {}
    pca_df_by_fold = {}

    if run_gwas:
        gwas_futures = []
        executor = None
        if Tools.blink_gwas_workers > 1:
            executor = ThreadPoolExecutor(max_workers=Tools.blink_gwas_workers)
            print(f"Running BLINK in parallel with {Tools.blink_gwas_workers} workers")
        for fold_index, (train_index, _) in enumerate(splits, start=1):
            print(f"Fold {fold_index}/{Tools.num_splits}: Preparing GWAS inputs...")
            train_genotype = genotype_df.iloc[train_index]
            train_genotype_numeric = Tools.convert_genotypes_to_numeric(train_genotype)
            principal_df = None

            if Tools.use_pcs and len(Tools.n_pcs_list) > 0:
                principal_df, explained_variance_ratio = Tools.perform_pca(
                    train_genotype_numeric,
                    fold_index,
                    n_components=max(len(Tools.n_pcs_list), 15),
                )
                explained_variances_by_fold[fold_index] = explained_variance_ratio
                pca_df_by_fold[fold_index] = principal_df

            file_fold_index = f"{fold_index}"
            for pc in Tools.n_pcs_list:
                print(f"Include {pc} PCs as covariant")
                file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
                Tools.blink_file_generator(train_genotype_numeric, file_fold_npc_index, principal_df, pc)
                if executor:
                    gwas_futures.append(executor.submit(Tools.run_gwas, file_fold_npc_index))
                else:
                    Tools.run_gwas(file_fold_npc_index)

        if executor:
            for future in as_completed(gwas_futures):
                future.result()
            executor.shutdown()

        if Tools.use_pcs and len(Tools.n_pcs_list) > 0 and pca_df_by_fold:
            Tools.pca_plotter(pca_df_by_fold)
            Tools.combined_scree_plot(explained_variances_by_fold)

    models_list = []
    snps_list = []
    if run_model:
        for fold_index, (train_index, test_index) in enumerate(splits, start=1):
            print(f"Fold {fold_index}/{Tools.num_splits}: Preparing model data...")
            train_genotype = genotype_df.iloc[train_index]
            test_genotype = genotype_df.iloc[test_index]
            file_fold_index = f"{fold_index}"
            for pc in Tools.n_pcs_list:
                file_fold_npc_index = f"{file_fold_index}_{pc}PCs"
                for n_snp in Tools.top_n_snps_list:
                    file_fold_npc_nsnp_index = f"{file_fold_npc_index}_{n_snp}SNPs"
                    current_dir = os.path.join(Tools.model_output_dir, file_fold_npc_nsnp_index)
                    os.makedirs(current_dir, exist_ok=True)
                    selected_snps = Tools.select_markers_and_prepare_data_for_model_construction(
                        file_fold_npc_index,
                        file_fold_npc_nsnp_index,
                        test_genotype,
                        top_n_snps=n_snp,
                    )
                    snps_list.append([file_fold_index, pc, n_snp, selected_snps])
                    # Prepare features
                    train_csv_path = os.path.join(current_dir, f'train_selected_snps_{file_fold_npc_nsnp_index}.csv')
                    test_csv_path = os.path.join(current_dir, f'test_selected_snps_{file_fold_npc_nsnp_index}.csv')

                    train_data = pd.read_csv(train_csv_path, index_col='taxa')
                    test_data = pd.read_csv(test_csv_path, index_col='taxa')

                    # Separate features and labels
                    train_y = train_data['Phenotype'].values.astype(int)
                    train_X = train_data.drop(columns=['Phenotype']).values

                    test_y = test_data['Phenotype'].values.astype(int)
                    test_X = test_data.drop(columns=['Phenotype']).values

                    feature_list = train_data.drop(columns=['Phenotype']).columns.tolist()

                    if Tools.oversampling:
                        smt = SMOTE()
                        train_X, train_y = smt.fit_resample(train_X, train_y)

                    for model_name in models:
                        models_list.append(
                            (model_name, train_X, train_y, test_X, test_y, file_fold_index, pc, n_snp, feature_list)
                        )
    
    if run_model and models_list:
        # Parallelize the machine learning process in parallel
        parallel_process(
                models_list,
                Tools.run_machine_learning_and_save_results,
                max_workers=Tools.num_threads
        )

    if run_model:
        # Concatenate results
        Tools.cat_result()
        # Collect SNPs from all folds
        snps_list_columns = ['file_fold_index', 'pc', 'n_snp', 'selected_snps']
        snps_list_df = pd.DataFrame(snps_list, columns=snps_list_columns)

        # Identify SNPs that appear more than once
        snps_list_file = os.path.join(Tools.output_dir, f'snps_list.txt')
        snps_list_df.to_csv(snps_list_file, index=False)

    print("Finsihed")
  
if __name__ == "__main__":
    run_snp_prediction_model()
