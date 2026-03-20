import pandas as pd
from typing import Optional, Dict
##dnnTL__FX_M01_CF_M04_MTLP_TL_ChengMJ
import sys
import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
DevicePath = "D:/Project/0.MutClone/Function"
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/BaseEstimator/2.1_dnnTL_multi_source_labeled_target')
from MultiSource_TargetLabel_dnnTL import * #import the neural-network classifier function for classifying source and target domains
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code')
from Compute_Transferability_Between_Source_and_Target_Domains_DTE_GPU import *

#Bayesian hyperparameter optimization
import optuna
from optuna.samplers import TPESampler, RandomSampler, CmaEsSampler, GridSampler
from optuna.pruners import MedianPruner, SuccessiveHalvingPruner, HyperbandPruner #pruners
import copy

#random
import random


#file I/O
import pickle
import os
import json
import torch
'''
@Point Python implementation of this framework
source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[统一数据类型]]/分类模型/[统一标准]分类模型-1.1.2.R.R")

 '''
class ClassifyModel:
    def __init__(self, 
                 data: pd.DataFrame=None, 
                 sample_label: pd.DataFrame = None,
                 classifier: str = None,
                 cross_validation_models: Optional['CreateClassifyModel'] = None,
                 CV_PP: dict = {'score':None,
                               'label':None,
                               'AUROC':None,
                               'AUPR':None},
                 final_model: Optional['MultiSource_TargetLabel_dnnTL'] = None,
                 #contents related to the standard data center
                 data_path_rds: str = None,
                 data_name: str = "APC",
                 way: str = "mut-exprs",
                 background: bool = False ) :
        super(ClassifyModel, self).__init__()
        #initialize attributes
        self.data = data
        self.sample_label = sample_label
        self.cross_validation_models = cross_validation_models
        self.CV_PP = CV_PP
        
        #initialize classifier
        self.classifier = classifier
        self.init_cm_Estimator()
        
        #check data; if missing, fill based on data_path_rds
        self.data_path_rds = data_path_rds
        self.data_name = data_name
        self.way = way
        self.background = background
        self.extract_data_label()
    
    def extract_data_label(self, way=None, background=None):
        #if parameters are not specified in the method, use the corresponding class attributes
        way = way or self.way
        background = background or self.background
        
        # data is directly provided in the model    
        if self.data is not None and self.sample_label is not None:
            print("数据直接被提供\n")
        else:
            if self.data_path_rds is not None and self.data_name is not None:
                if way == "mut-exprs":  # use standard data storage format
                    # read data from the RDS file
                    with open(self.data_path_rds, 'rb') as file:
                        t_data = pickle.load(file)
                        all_data = t_data['all.data']
                        all_label = t_data['all.label']
                        temp = all_label[self.data_name]
                    
                    if background:
                        self.data = all_data.loc[temp['background'], temp['label'].keys()]
                    else:
                        self.data = all_data.loc[:, temp['label'].keys()]
                        
                    self.sample_label = pd.DataFrame.from_dict(temp['label'], orient='index')
                    print("从数据地址中完成数据导入\n")   
            else:
                raise ValueError("没有取到该模型相关数据: 数据地址或数据全局变量有问题...")
    
    
    #(0) Initialize classifier
    def init_cm_Estimator(self, **kwargs):
        if self.classifier == "MultiSource_TargetLabel_dnnTL":
            self.final_model = MultiSource_TargetLabel_dnnTL(**kwargs)
    '''
    (1) Call GridSearchCV to search the hyperparameter grid'''
    # param_grid = {
    #     'hiddenFX_dim': [128, 256, 512],
    #     'hiddenCF_dim': [128, 256, 512],
    #     'dropout_rate': [0.3, 0.5, 0.7],
    #     'batch_size': [16, 32, 64],
    #     'epoch': [50, 100],
    #     'lr': [1e-4, 1e-3, 1e-2],
    #     'lam1': [0.1, 0.2, 0.3],
    #     'n_source':[1,2,3]}
    def HyperparasGridSearchCV(self, param_grid=None, cv=5, threads=10, model=None, GPU_id=None):
        #set the search grid and related parameters
        self.final_model.GPU_id = GPU_id
        grid_search = GridSearchCV(estimator=self.final_model, param_grid=param_grid, scoring='roc_auc', cv=cv, n_jobs=threads) #, error_score='raise'
        grid_search.fit(X=self.data.T, y=self.sample_label, model=model) #start searching
        #store the best search result
        self.OptimalHyperparameters = grid_search.best_params_
        return(grid_search.best_params_)
    
    def generate_valid_cv_splits(self, cv=5, max_retry=1000):
        """
        ✅ Automatically generate stratified k-fold splits without single-class folds
        """
        labels = self.sample_label
        for attempt in range(1, max_retry + 1):
            skf = StratifiedKFold(n_splits=cv,shuffle=True,random_state=random.randint(0, 999999))
            folds = []
            bad_split = False
            for train_idx, test_idx in skf.split(self.data.T, labels):
                y_train = labels.iloc[train_idx]
                y_test = labels.iloc[test_idx]
                if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
                    bad_split = True
                    break
                folds.append((train_idx, test_idx))
            if not bad_split:
                print(f"✅ CV 划分成功（尝试 {attempt} 次）")
                return folds
        raise RuntimeError(
            f"❌ 尝试 {max_retry} 次后仍无法获得合法 CV 划分，请检查数据类别分布"
        )


    # -------------------- Added Optuna hyperparameter optimization --------------------
    def HyperparasOptunaCV(
        self,
        param_space: Dict[str, dict],
        cv: int = 5, 
        n_trials: int = 30, #number of sampling trials
        threads: int = 1,
        model=None,
        GPU_id=None,
        sampler_name: str = "tpe",
        pruner_name: str = "median",   # added: pruning strategy
        random_state: int = 8766
    ):
        """
        Use Optuna for hyperparameter optimization, accepting two input forms:
        1️) Range-based input: {"type": "int", "low": 128, "high": 512, "step": 128}
        2️) List-based input: {"type": "int", "values": [128, 256, 512]}
        sampler_name: "tpe" | "random" | "cmaes" | "grid"

        #input examples
        #example1
        param_space = {
                "hiddenFX_dim": {"type": "int", "low": 128, "high": 512, "step": 128},
                "dropout_rate": {"type": "float", "low": 0.3, "high": 0.7},
            }

        #example2
        param_space = {
                "hiddenFX_dim": {"type": "int", "values": [128, 256, 512]},
                "batch_size": {"type": "int", "values": [16, 32, 64]},
            }

        #example2
        param_space = {
                "hiddenFX_dim": {"type": "int", "low": 128, "high": 512, "step": 128},
                "batch_size": {"type": "int", "values": [16, 32, 64]},
            }
        """
        # select different samplers
        if sampler_name.lower() == "random":
            sampler = RandomSampler(seed=random_state)
        elif sampler_name.lower() == "cmaes":
            sampler = CmaEsSampler(seed=random_state)
        elif sampler_name.lower() == "grid":
            sampler = GridSampler(search_space=param_space)
        else:
            sampler = TPESampler(seed=random_state)

        # choose pruner
        if pruner_name.lower() == "median":
            pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=2)
        elif pruner_name.lower() == "halving":
            pruner = SuccessiveHalvingPruner()
        elif pruner_name.lower() == "hyperband":
            pruner = HyperbandPruner()
        else:
            pruner = optuna.pruners.NopPruner()  # no pruning
    
        study = optuna.create_study(direction="maximize", sampler=sampler, pruner=pruner)
        
        
        # ------------------------- 2) Fix the split once (shared by all trials) -------------------------
        counts = np.unique(self.sample_label, return_counts=True)[1]
        if np.min(counts) < cv:
            raise ValueError(f"类别数量不足以进行 {cv}-折交叉验证")
        folds = self.generate_valid_cv_splits(cv=cv)
        

        # ------------------------- 3) Record details of the best trial -------------------------
        best_trial_data = {"score": None, "label": None, "AUROC": 0, "AUPR": 0}
        # define the objective function
        def objective(trial):
            params = {}
            for key, p in param_space.items():
                ptype = p.get("type", "float").lower()
    
                # ---- support list-based input ----
                if "values" in p:
                    params[key] = trial.suggest_categorical(key, p["values"])
    
                # ---- support range-based input ----
                elif ptype == "int":
                    params[key] = trial.suggest_int(key, p["low"], p["high"], step=p.get("step", 1))
                elif ptype == "float":
                    params[key] = trial.suggest_float(key, p["low"], p["high"], step=p.get("step", None)) #if step is None, it is continuous
                elif ptype == "loguniform":
                    params[key] = trial.suggest_float(key, p["low"], p["high"], log=True)
                elif ptype == "categorical":
                    params[key] = trial.suggest_categorical(key, p["choices"])
                else:
                    raise ValueError(f"未知参数类型或格式: {p}")
    
            # update parameters into the model instance
            #self.final_model.__dict__.update(params)
            #self.final_model.GPU_id = GPU_id
    
            # start cross-validation using the current split #compute mean AUROC by cross-validation
            auroc_scores, aupr_scores = [], []
            fold_scores, fold_labels = {}, {}
            for fold_idx, (train_idx, test_idx) in enumerate(folds):
                # update parameters into the model instance
                model_copy = copy.deepcopy(self.final_model) #a copy is required so that each hyperparameter setting starts from the same initial point
                model_copy.GPU_id = GPU_id
                model_copy.__dict__.update(params)

                X_train, X_test = self.data.iloc[:, train_idx], self.data.iloc[:, test_idx]
                y_train, y_test = self.sample_label.iloc[train_idx], self.sample_label.iloc[test_idx]

                model_i = model_copy.fit(X=X_train.T, y=y_train, model=model)
                preds = model_i.Predict(X=X_test.T)

                # convert format for unified storage
                test_sample_names = self.data.columns[test_idx]
                predictions_df = pd.DataFrame(preds.values, index=test_sample_names, columns=["score"])
                label_df = pd.DataFrame(y_test.values, index=test_sample_names, columns=["label"])
                fold_scores[fold_idx] = predictions_df
                fold_labels[fold_idx] = label_df

                # compute metrics
                auroc = roc_auc_score(y_test, preds)
                aupr = average_precision_score(y_test, preds)
                auroc_scores.append(auroc)
                aupr_scores.append(aupr)

                trial.report(np.mean(auroc_scores), step=fold_idx)
                if trial.should_prune():
                    raise optuna.TrialPruned()

            mean_auroc = np.mean(auroc_scores)
            mean_aupr = np.mean(aupr_scores)

            # if this trial is better than the current best, record cross-validation information
            if mean_auroc > best_trial_data["AUROC"]:
                best_trial_data.update({
                    "score": fold_scores,
                    "label": fold_labels,
                    "AUROC": mean_auroc,
                    "AUPR": mean_aupr
                })

            return mean_auroc if len(auroc_scores) > 0 else 0.0

        # execute search
        study.optimize(objective, n_trials=n_trials, n_jobs=threads, show_progress_bar=True, catch=(RuntimeError,))

        # record the cross-validation results of the best trial
        self.CV_PP = best_trial_data     
    
        # save the best parameters and score
        self.OptimalHyperparameters = study.best_params
        self.final_model.__dict__.update(self.OptimalHyperparameters)
        self.best_optuna_score = study.best_value
        print(f"Optuna搜索完成 | 最优AUROC: {study.best_value:.4f}")
        print(f"最优参数: {study.best_params}")
    
        return study.best_params, study.best_value


    #(2) Perform cross-validation evaluation based on the best hyperparameters and fill CV_PP
    def CV_eval(self, cv=5, GPU_id=None, **kwargs):
        #skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=8766)
        auroc_scores = []
        aupr_scores = []
        fold_scores = {}
        fold_labels = {}
        self.final_model.GPU_id = GPU_id

        # obtain valid splits
        counts = np.unique(self.sample_label, return_counts=True)[1]
        if np.min(counts) < cv:
            raise ValueError(f"类别数量不足以进行 {cv}-折交叉验证")
        folds = self.generate_valid_cv_splits(cv=cv)
        
        for fold, (train_index, test_index) in enumerate(folds):
            # update the best hyperparameters into the initial model instance
            model_copy = copy.deepcopy(self.final_model) #a copy is required so that each hyperparameter trail starts from the same initial model 
            model_copy.GPU_id = GPU_id
            model_copy.__dict__.update(self.OptimalHyperparameters)
            #data
            X_train, X_test = self.data.iloc[:, train_index], self.data.iloc[:, test_index]
            y_train, y_test = self.sample_label.iloc[train_index], self.sample_label.iloc[test_index]
            #training
            model = model_copy.fit(X=X_train.T, y=y_train, **kwargs) #train the model on the training set with the best hyperparameters
            predictions = model.Predict(X=X_test.T) 

            #unify data format
            test_sample_names = self.data.columns[test_index]
            predictions_df = pd.DataFrame(predictions.values, index=test_sample_names, columns=["score"])
            label_df = pd.DataFrame(y_test.values, index=test_sample_names, columns=["label"])
            
            fold_scores[fold] = predictions_df
            fold_labels[fold] = label_df
            
            auroc_scores.append(roc_auc_score(y_test, predictions))
            aupr_scores.append(average_precision_score(y_test, predictions))
        
        self.CV_PP['score'] = fold_scores
        self.CV_PP['label'] = fold_labels
        self.CV_PP['AUROC'] = np.mean(auroc_scores)
        self.CV_PP['AUPR'] = np.mean(aupr_scores)
    
    #(3) Train the model based on the class attribute OptimalHyperparameters
    def Train_best_model(self, GPU_id=None, **kwargs):
        self.final_model.GPU_id = GPU_id
        opth = self.OptimalHyperparameters
        self.final_model.__dict__.update(opth)  #fill the model with the best hyperparameters
        best_model = self.final_model.fit(X=self.data.T, y=self.sample_label, **kwargs)  #input the dictionary into the best model
        self.final_model = best_model  #store the output in the final_model attribute of the current ClassifyModel
        return(best_model)
    
    
    #(4) Add a function for calculating the optimal threshold
    def compute_best_threshold(self, pos_label, best_method="youden"):
        all_scores = []
        all_labels = []
        
        # merge score and label from all folds
        for fold, score_df in self.CV_PP['score'].items():
            all_scores.extend(score_df['score'].values)
            all_labels.extend(self.CV_PP['label'][fold]['label'].values)
        
        # convert to numpy arrays
        all_scores = np.array(all_scores)
        all_labels = np.array(all_labels)
        
        # call RocGetThreshold to calculate the optimal threshold
        best_threshold = RocGetThreshold(score=all_scores, label=all_labels, pos=pos_label, best_method=best_method)
        self.optimal_threshold = best_threshold
        
        return best_threshold
    
    #(5) Final-model prediction
    def final_model_Predict(self, predict_data, normalizeUsingCM=False, imputeZero=False, GPU_id=None):
        #predict_data rows are genes and columns are samples
        self.extract_data_label()           #extract data and labels into the classification object
        
        #ensure that predict_data and training data have the same genes
        common_genes = self.data.index.intersection(predict_data.index)
        if common_genes.empty:
            raise ValueError("No common genes between training data and predict_data.")
        
        train_data_common = self.data.loc[common_genes]
        predict_data_common = predict_data.loc[common_genes]
        
        #if quantile normalization is used
        if normalizeUsingCM:
            predict_data_common = normalize_quantile_cross_platform(ref_matx=train_data_common, query_matx=predict_data_common)
        
        #if missing values are imputed with zeros
        if imputeZero:
            features = self.final_model.model.features
            # Create a matrix of zeros with rows as features and columns as predict_data columns
            imputed_data = pd.DataFrame(0.0, index=features, columns=predict_data.columns)
            # Fill with existing data for common features
            common_features = predict_data_common.index.intersection(features)
            imputed_data.loc[common_features] = predict_data_common.loc[common_features]
            predict_data = imputed_data
                    
        y = self.final_model.Predict(X=predict_data.T, GPU_id=GPU_id)
        res = y
        # Check if self.optimal_threshold exists and apply it
        if hasattr(self, 'optimal_threshold') and self.optimal_threshold is not None:
            labels = (y >= self.optimal_threshold).astype(int)
            res = pd.DataFrame({'Prediction': y.values.flatten().tolist(), 'label': labels.values.flatten().tolist()}, index=y.index.astype(str))                 
            
        return(res)
    
    
    #' @description: Export final model weights, best hyperparameters, and the basic/summary information of ClassifyModel for subsequent reproduction and inspection
    #' @param export_dir [str] export directory, e.g., "export_model"
    #' @return [None]
    #  @code logic
    #  @(1) save the weights of final_model.model (.pt)
    #  @(2) save the best hyperparameters (.json)
    #  @(3) save a ClassifyModel information dictionary containing only basic fields (including CV_PP, .pkl), which can be loaded without depending on the ClassifyModel class definition
    #  @(4) save summary information of ClassifyModel without CV_PP (.json)
    def export(self, export_dir="export_model"):
        """
        Export all essential model information:
        -------------------------------------------------
        (1) weights of final_model.model → final_model_state.pt
        (2) best model hyperparameters → hyperparameters.json
        (3) basic ClassifyModel information dictionary (including CV_PP) → ClassifyModel_basic.pkl
        (4) ClassifyModel summary information (excluding CV_PP) → ClassifyModel_summary.json
        -------------------------------------------------
        """
        # create export directory (automatically if it does not exist)
        os.makedirs(export_dir, exist_ok=True)

        # =======================================================
        # 1) Save model weights (final_model.model.state_dict)
        #    Compatible with both single-GPU and DataParallel multi-GPU models
        # =======================================================
        model_state_path = os.path.join(export_dir, "final_model_state.pt")
        if hasattr(self, "final_model") and hasattr(self.final_model, "model") and self.final_model.model is not None:
            model_to_save = self.final_model.model
            # if DataParallel is used, take the internal module
            if isinstance(model_to_save, torch.nn.DataParallel):
                model_to_save = model_to_save.module
            #! uniformly map model parameters to CPU before saving to ensure safe loading in CPU-only environments
            torch.save(model_to_save.cpu().state_dict(), model_state_path)
            print(f"✔ 模型权重已保存：{model_state_path}")
        else:
            print("⚠ 无 final_model.model，无法保存模型权重")

        # =======================================================
        # 2) Save hyperparameters (json)
        # =======================================================
        hyper_params = getattr(self, "OptimalHyperparameters", None)
        if hyper_params is None:
            hyper_params = {}
            print("⚠ 当前对象中不存在 OptimalHyperparameters, 导出空超参数字典")

        hyper_param_path = os.path.join(export_dir, "hyperparameters.json")
        with open(hyper_param_path, "w", encoding="utf-8") as f:
            json.dump(hyper_params, f, indent=4, ensure_ascii=False)
        print(f"✔ 超参数已保存：{hyper_param_path}")

        # =======================================================
        # 3) Save basic ClassifyModel information (including CV_PP) → pkl
        #    Keep only the basic fields needed for downstream analysis, without directly saving class instances such as ClassifyModel / final_model
        # =======================================================
        # 3.1) Automatically keep all fields except final_model
        exclude_keys = {"final_model"}
        basic_info = {
            k: v
            for k, v in self.__dict__.items()
            if k not in exclude_keys
        }
        # 3.2) Add key information needed from final_model
        basic_info.update({
            "features": self.final_model.model.features,
            "init_hyperparameters": self.final_model.model.hyperparameters
        })

        full_pkl_path = os.path.join(export_dir, "ClassifyModel_basic.pkl")
        with open(full_pkl_path, "wb") as f:
            pickle.dump(basic_info, f)
        print(f"✔ ClassifyModel 基础信息（含 CV_PP）已保存：{full_pkl_path}")

        # =======================================================
        # 4) Output ClassifyModel summary information (json, excluding CV_PP)
        # =======================================================
        cv_pp = getattr(self, "CV_PP", {}) or {}
        summary = {
            "classifier": getattr(self, "classifier", None),
            "data_name": getattr(self, "data_name", None),
            "way": getattr(self, "way", None),
            "background": getattr(self, "background", None),
            "OptimalHyperparameters": hyper_params,
            "optimal_threshold": getattr(self, "optimal_threshold", None),
            "CV_AUROC": cv_pp.get("AUROC", None) if isinstance(cv_pp, dict) else None,
            "CV_AUPR": cv_pp.get("AUPR", None) if isinstance(cv_pp, dict) else None,
        }

        summary_path = os.path.join(export_dir, "ClassifyModel_summary.json")
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=4, ensure_ascii=False)
        print(f"✔ ClassifyModel 摘要已保存：{summary_path}")

        print("\n🎉 模型已全部导出完毕！")



'''
@Subfunction_estimate the optimal threshold based on cross-validation '''
# Define RocGetThreshold function
def RocGetThreshold(score, label, pos=None, best_method="youden", r=1.0):
    label = np.array(label)
    if pos is not None:
        label = np.where(label == pos, 1, 0)
        
    fpr, tpr, thresholds = roc_curve(label, score)
    
    if best_method == "youden":
        youden_index = tpr + r * (1 - fpr) - 1
        best_threshold_index = np.argmax(youden_index)
    elif best_method == "closest.topleft":
        distances = np.sqrt((1 - tpr) ** 2 + fpr ** 2)
        best_threshold_index = np.argmin(distances)
    else:
        raise ValueError(f"未知的最佳方法: {best_method}")
    
    return thresholds[best_threshold_index]


'''
@Subfunction_quantile-normalize a new dataset using the training data as reference '''
#import numpy as np
#import pandas as pd
def normalize_quantile_cross_platform(ref_matx, query_matx):
    """
    Perform cross-platform quantile normalization on query matrix using the reference matrix distribution.
    Parameters:
        ref_matx (pd.DataFrame): Reference matrix [genes x samples].
        query_matx (pd.DataFrame): Query matrix [genes x samples] to normalize based on ref_matx distribution.

    Returns:
        numpy.ndarray: Quantile normalized query matrix.
    """
    # Point[1] Calculate the quantile target from the reference matrix
    sorted_ref = np.sort(ref_matx.values, axis=0)  # Sort each column of the reference matrix
    qn_target = np.mean(sorted_ref, axis=1)       # Compute the mean for each quantile across columns

    # Point[2] Apply quantile normalization to the query matrix
    sorted_query = np.sort(query_matx.values, axis=0)  # Sort each column of the query matrix
    ranks = np.argsort(np.argsort(query_matx.values, axis=0), axis=0)  # Get the ranks for each element in the query matrix
    qn_query = np.zeros_like(query_matx.values, dtype=np.float64)
    for col in range(query_matx.shape[1]):
        qn_query[:, col] = qn_target[ranks[:, col]]
        
    # Convert the result back to a DataFrame with the original index and column names
    qn_query_df = pd.DataFrame(qn_query, index=query_matx.index, columns=query_matx.columns)
    return qn_query_df

