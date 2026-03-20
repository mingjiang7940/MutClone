import os
import json
import pickle
import torch
import sys
DevicePath = "D:/Project/0.MutClone/Function"
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/BasicArchitecture')
'''
@Important_Do_not_change_the_function_filename_lightly_if_changed_previously_saved_trained_model_pkl_files_will_not_be_imported_correctly_into_the_python_environment'''
from FX001_MTLP001_MTLP002 import *

def load_ClassifyModel(export_dir):  #, ClassifyModel_class, model_class
    """
    Restore a complete and directly usable ClassifyModel object from the 2 exported files.
    "final_model_state.pt"
    "ClassifyModel_basic.pkl"
    
    Parameters
    -----
    export_dir : str
        Directory containing the exported files
    ClassifyModel_class : classification model class
        For example, ClassifyModel 
        
    model_class : the final trained neural-network framework that can directly receive the neural-network weight parameters [the class of cm.final_model.model]
        For example, FX001_MTLP001_MTLP002

    Returns
    -----
    cm : ClassifyModel object
        Restored state, can be used directly for predict()
        
    #Running example:
    #Important: do not put quotation marks around the two class-name parameters, because they directly refer to classes (just use them directly after import; if the file or class name changes, update the import accordingly)
    cm = load_ClassifyModel(export_dir="/WorkSpace/chengmingjiang/TmpData/TL_model/CRC_TCGA/2.迁移学习dnnTL/单个基因/TP53", 
                            ClassifyModel_class=ClassifyModel, model_class=FX001_MTLP001_MTLP002)   
    """
    #01) model weights
    model_state_path = os.path.join(export_dir, "final_model_state.pt")
    state_dict = torch.load(model_state_path, map_location="cpu")
    #02) basic information (independent of ClassifyModel)
    basic_path = os.path.join(export_dir, "ClassifyModel_basic.pkl")
    with open(basic_path, "rb") as f:
        basic_info = pickle.load(f)

    # ====================================================
    # 2)          Initialize class model_class from the basic model information basic_info
    # add_extra_fields: whether to store other attributes in basic_info into the initialized class in addition to the attributes valid for initialization
    # ====================================================
    import inspect
    def init_from_basic(model_class, basic_info, add_extra_fields=True):
        """
        Initialize model_class using basic_info
        - Only fields supported by __init__ are used
        - Optionally attach extra fields to the model object
        """
        # 1) First extract valid init fields
        #init_kwargs = filter_valid_hyperparams(basic_info, model_class)
        sig = inspect.signature(model_class.__init__)
        valid_keys = set(sig.parameters.keys()) - {"self"}
        init_kwargs = {k: v for k, v in basic_info.items() if k in valid_keys}
        
        # 2) Initialize the model (using only valid fields)
        model = model_class(**init_kwargs)
        # 3) Optionally add non-init fields to the model
        if add_extra_fields:
            for k, v in basic_info.items():
                if k not in init_kwargs:
                    setattr(model, k, v)
        return model

    if basic_info['classifier']=='MultiSource_TargetLabel_dnnTL':
        ClassifyModel_class = ClassifyModel
        cm = init_from_basic(ClassifyModel_class, basic_info, add_extra_fields=True)    #extra fields need to be added here because these are post-training parameters
        cm.final_model.__dict__.update(basic_info["OptimalHyperparameters"])

        # ====================================================
        # 3) Build final_model.model neural network and load weights
        # ====================================================
        # Initialize neural-network architecture
        model_class = FX001_MTLP001_MTLP002
        cm.final_model.model = init_from_basic(model_class, basic_info["init_hyperparameters"],add_extra_fields=False)    #do not add extra fields here; the other fields are already covered by the optimal hyperparameters and do not need to be repeatedly added
        cm.final_model.model.features = basic_info["features"]   #this must be added for use
        cm.final_model.model.load_state_dict(state_dict)         #load weights
        cm.final_model.model.hyperparameters = basic_info["init_hyperparameters"]
    
    return cm

