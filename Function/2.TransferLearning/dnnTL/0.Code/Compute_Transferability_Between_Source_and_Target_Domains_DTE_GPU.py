import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import minmax_scale
''' Python implementation of DTE-based source ranking
source("/pub5/xiaoyun/BioX/DataScience/机器学习/[深度学习]/Learning迁移学习/[1.源的选择]/DTE[可用于迁移学习选源]/计算源域和靶域的可迁移性DTE.R")

'''
def calculate_scatter_matrix(dat_matrix, lab, GPU_id=None):
    device = torch.device(f"cuda:{GPU_id}" if GPU_id is not None else "cpu")
    """
    Compute the between-class scatter matrix (on GPU)
    :param dat_matrix: data matrix
    :param lab: labels
    :param device: 'cpu' or 'cuda'
    :return: between-class scatter matrix
    """
    # Convert data to PyTorch tensors and move them to GPU
    #dat_matrix = torch.tensor(dat_matrix, dtype=torch.float32, device=device)
    dat_matrix = dat_matrix.clone().detach().to(device).float()
    lab = torch.tensor(lab, dtype=torch.float32, device=device)
    
    # Compute the centroid of each class
    unique_labels = torch.unique(lab)
    mean_each_class = torch.stack([dat_matrix[:, lab == l].mean(dim=1) for l in unique_labels], dim=1)

    # Compute the global centroid of the data
    mean_total = dat_matrix.mean(dim=1, keepdim=True)
    
    # Compute the number of samples in each class
    samp_num_each_class = torch.tensor([torch.sum(lab == l).item() for l in unique_labels], device=device)
    
    # Compute the deviation matrix
    tpm_deviation_matrix = mean_each_class - mean_total

    # For each class, multiply the deviation vector by its transpose to compute the class-specific Sb matrix
    Sb_matrix_each_class = [torch.ger(tpm_deviation_matrix[:, i], tpm_deviation_matrix[:, i]) 
                            for i in range(tpm_deviation_matrix.shape[1])]
    
    # Multiply by the sample size of each class
    Sb_matrix_each_class = [Sb * samp_num_each_class[i] for i, Sb in enumerate(Sb_matrix_each_class)]

    # Sum the between-class scatter matrices of all classes
    Sbetween = sum(Sb_matrix_each_class)

    return Sbetween

def calculate_DIS_DIF(source_data, source_lab, target_data, GPU_id=None):
    device = torch.device(f"cuda:{GPU_id}" if GPU_id is not None else "cpu")
    """
    Compute the DIS and DIF metrics (on GPU)
    :param source_data: source-domain data
    :param source_lab: source-domain labels
    :param target_data: target-domain data
    :param device: 'cpu' or 'cuda'
    :return: DIS and DIF values
    """
    # Trim data
    commgene = source_data.index.intersection(target_data.index)
    source_data = source_data.loc[commgene]
    target_data = target_data.loc[commgene]

    # Move data to GPU
    source_data_tensor = torch.tensor(source_data.values, dtype=torch.float32, device=device)
    target_data_tensor = torch.tensor(target_data.values, dtype=torch.float32, device=device)

    # Compute the between-class scatter matrix of the source domain
    #print("Compute between-class scatter matrix")
    Sbs = calculate_scatter_matrix(source_data_tensor, source_lab, GPU_id=GPU_id)

    # Compute the between-domain scatter matrix of source and target
    #print("Compute between-domain scatter matrix of source and target")
    combined_data = torch.cat([source_data_tensor, target_data_tensor], dim=1)
    # Convert labels to numeric type
    combined_labels = np.array([0] * source_data.shape[1] + [1] * target_data.shape[1])
    Sbst = calculate_scatter_matrix(combined_data, combined_labels, GPU_id=GPU_id)

    # Compute Domain Transferability Estimation
    print("Compute DIS and DIF in DTE")
    # DIS = torch.norm(Sbs, p=1).item()
    # DIF = torch.norm(Sbst, p=1).item()

    # We first need to compute the sum of absolute values for each column
    abs_sum_columns_Sbs = torch.sum(torch.abs(Sbs), dim=0)
    abs_sum_columns_Sbst = torch.sum(torch.abs(Sbst), dim=0)
    # Then take the maximum of these sums
    DIS = torch.max(abs_sum_columns_Sbs).item()
    DIF = torch.max(abs_sum_columns_Sbst).item()
    # DIS = torch.max(torch.sum(torch.abs(Sbs), dim=0)).item()
    # DIF = torch.max(torch.sum(torch.abs(Sbst), dim=0)).item()

    return DIS, DIF

def sort_source_by_DTE(source_data_list, source_lab_list, target_data, GPU_id=None):
    """
    Rank source domains according to the DTE metric (on GPU)
    :param source_data_list: dictionary of source-domain data
    :param source_lab_list: dictionary of source-domain labels
    :param target_data: target-domain data
    :param GPU_id: 'cpu' or 'cuda'
    :return: sorted DTE dictionary
    """
    DIS_DIF_list = {name: calculate_DIS_DIF(source_data_list[name], source_lab_list[name], target_data, GPU_id=GPU_id) 
                    for name in source_data_list.keys()}

    DIS_ls = np.array([DIS_DIF[0] for DIS_DIF in DIS_DIF_list.values()])
    DIF_ls = np.array([DIS_DIF[1] for DIS_DIF in DIS_DIF_list.values()])

    # Normalize
    DIS_ls = minmax_scale(DIS_ls, feature_range=(0, 1))
    DIF_ls = 1 - minmax_scale(DIF_ls, feature_range=(0, 1))

    DTE = DIS_ls * DIF_ls
    DTE_dict = {name: DTE[i] for i, name in enumerate(DIS_DIF_list.keys())}

    DTE_sorted = dict(sorted(DTE_dict.items(), key=lambda item: item[1], reverse=True))

    return DTE_sorted



