import os
import pickle
import pandas as pd
import sys
import time
DevicePath = "D:/Project/0.MutClone/Function"
sys.path.append(DevicePath + '/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *


'''
@构建dnnTL突变预测模型
'''
ID_xteam = "CRC_TCGA"
out_dir = os.path.join(os.path.dirname(DevicePath), "0.Data/Results", ID_xteam, "3.dnnTL")
os.makedirs(out_dir, exist_ok=True)

way = "mut-exprs"
background = False
method_select_source = "DTE"

nthreads = 2
GPU_id = 0

pkl_file_path = os.path.join(os.path.dirname(DevicePath), "0.Data","CancerDriverGene.pkl")
with open(pkl_file_path, 'rb') as pkl_file:
    pancancer_drive = pickle.load(pkl_file)

with open( os.path.join(os.path.dirname(DevicePath), "0.Data/TCGA_MutClone", ID_xteam, 'mut.pancancer.exprs.data.indel.pkl'), 'rb') as file:
    all_label = pickle.load(file)['all.label']

mutgenes = list(set(all_label.keys()) & set(pancancer_drive))
#mutgenes = ["BOD1L1"]
gene_dir = os.path.join(out_dir, "Models")
os.makedirs(gene_dir, exist_ok=True)

ID_xteams = [
    "ACC_TCGA", "BLCA_TCGA", "BRCA_TCGA", "CESC_TCGA", "CHOL_TCGA", "CRC_TCGA",
    "DLBC_TCGA", "ESCA_TCGA", "GBM_TCGA", "HNSC_TCGA", "KICH_TCGA", "KIRC_TCGA",
    "KIRP_TCGA", "LGG_TCGA", "LIHC_TCGA", "LUAD_TCGA", "LUSC_TCGA", "MESO_TCGA",
    "OV_TCGA", "PAAD_TCGA", "PCPG_TCGA", "PRAD_TCGA", "SARC_TCGA", "SKCM_TCGA",
    "STAD_TCGA", "TGCT_TCGA", "THCA_TCGA", "THYM_TCGA", "UCEC_TCGA", "UCS_TCGA",
    "UVM_TCGA"
]
ID_xteams = list(set(ID_xteams) - {ID_xteam})

base_path = os.path.join(os.path.dirname(DevicePath), "0.Data/TCGA_MutClone")
aim_path = 'mut.pancancer.exprs.data.indel.pkl'

all_model = {}

if len(mutgenes) == 0:
    print('没有突变基因满足要求\n')
else:
    start_time = time.time()
    for mutgene in mutgenes:

        gene_file_path = os.path.join(gene_dir, f"{mutgene}.pkl")
        if os.path.exists(gene_file_path):
            print(f"跳过已存在的基因: {mutgene}")
            continue

        # (0) identify transferable source domains
        valid_ID_xteams = []
        for tmp_ID_xteam in ID_xteams:
            file_path = os.path.join(base_path, tmp_ID_xteam, aim_path)
            if os.path.exists(file_path):
                try:
                    with open(file_path, 'rb') as file:
                        all_label = pickle.load(file)['all.label']
                        if len(all_label) != 0 and mutgene in list(all_label.keys()):
                            valid_ID_xteams.append(tmp_ID_xteam)
                except (OSError, IOError, pickle.UnpicklingError) as e:
                    print(f"Error reading file for {tmp_ID_xteam}: {e}")

        if len(valid_ID_xteams) == 0:
            continue

        # (1) build the source models
        s_models = {}
        for tmp_ID_xteam in valid_ID_xteams:
            s_models[tmp_ID_xteam] = ClassifyModel(
                data_path_rds=os.path.join(base_path, tmp_ID_xteam, aim_path),
                data_name=mutgene, classifier="MultiSource_TargetLabel_dnnTL",
                way=way, background=background
            )
            s_models[tmp_ID_xteam].init_cm_Estimator(method_sort_sources=method_select_source)

        # (2) build the target model
        t_model = ClassifyModel(
            data_path_rds=os.path.join(base_path, ID_xteam, aim_path),
            data_name=mutgene, classifier="MultiSource_TargetLabel_dnnTL",
            way=way, background=background
        )
        t_model.init_cm_Estimator(method_sort_sources=method_select_source)

        # (3) hyperparameter optimization in cross-validation
        param_space = {
            "hiddenFX_dim": {"type": "int", "values": [256, 512]},
            "hiddenCF_dim": {"type": "int", "values": [128, 256]},
            "dropout_rate": {"type": "float", "values": [0.3, 0.5]},
            "batch_size": {"type": "int", "values": [16, 32]},
            "epoch": {"type": "int", "values": [50, 100]},
            "lr": {"type": "loguniform", "low": 1e-6, "high": 1e-2},
            "gamma": {"type": "float","low": 1, "high": 5},  
            "n_source": {"type": "int", "values": list(range(1, len(s_models) + 1))}
        }


        best_params, best_score = t_model.HyperparasOptunaCV(
            param_space=param_space,
            model=s_models,
            n_trials=30,
            sampler_name="tpe",
            cv=5,
            GPU_id=GPU_id,
            threads=nthreads
        )


        # (4) Train_best_model
        t_model.Train_best_model(model=s_models, GPU_id=GPU_id)

        # 保存模型
        t_model.data = None
        t_model.sample_label = None
        t_model.final_model.GPU_id = None
        base_dir2 = os.path.join(out_dir, "单个基因")
        t_model.export(export_dir=os.path.join(base_dir2, f"{gene}"))
        print(f"✅ 完成基因: {mutgene}\n")

print("🎉 所有基因建模完成！")
