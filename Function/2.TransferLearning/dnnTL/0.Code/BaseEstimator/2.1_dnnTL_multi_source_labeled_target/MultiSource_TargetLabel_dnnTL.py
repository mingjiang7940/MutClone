import torch
import numpy as np
from numpy import random
import pandas as pd
import sklearn.preprocessing as sk
from torch.utils.data import DataLoader, TensorDataset
import sys
from torch.utils.data.sampler import WeightedRandomSampler
import itertools
from itertools import cycle
from typing import Optional
from collections import OrderedDict
from _ast import For
#import pdb  #used for breakpoint debugging; after setting breakpoints, variables during function execution can be inspected

DevicePath = "D:/Project/0.MutClone/Function"
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/BasicArchitecture')
'''
@Important_Do_not_change_the_function_filename_lightly_if_changed_previously_saved_trained_model_pkl_files_will_not_be_imported_correctly_into_the_python_environment'''
from FX001_MTLP001_MTLP002 import *

sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code')
from Compute_Transferability_Between_Source_and_Target_Domains_DTE_GPU import *
#torch.set_num_threads(1)  #number of CPU threads used by pytorch; this must be set, otherwise our server will default to using 52 threads

##
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score


'''
@todo: Training and prediction of FX_M01_CF_M04_MTLP_TL_ChengMJ [inherits from BaseEstimator]

This class has 3 methods:
                  fit            model training method
                  Predict        model prediction method
                  predict_proba  model prediction method [same as Predict, but its output must be in array format, cannot carry sample names, and is only used for the hyperparameter search interface]
'''
class MultiSource_TargetLabel_dnnTL(BaseEstimator, ClassifierMixin):
    def __init__(self, 
                #[neural network initialization] hyperparameters
                hiddenFX_dim: int = 512,  #hidden-layer dimension of FX
                hiddenCF_dim: int = 256,  #hidden-layer dimension of CF
                output_dim: int = 1,
                dropout_rate: float = 0.5,
                
                #[model training] hyperparameters
                batch_size: int=32,
                epoch: int=100,
                lr: float=1e-3,        #learning rate
                n_source: int=1,       #number of sources
                method_sort_sources: str="DTE", #if None, do not rank sources
                gamma: float = 1.5,
                GPU_id: Optional[int] = None
                
        ):
        #[neural network initialization] hyperparameters
        self.hiddenFX_dim = hiddenFX_dim
        self.hiddenCF_dim = hiddenCF_dim
        self.output_dim = output_dim
        self.dropout_rate = dropout_rate
        
        #[model training] hyperparameters
        self.batch_size = batch_size
        self.epoch = epoch
        self.lr = lr
        self.n_source = n_source
        self.method_sort_sources = method_sort_sources

        #general class parameters for binary classification
        self.classes_ = np.array([0, 1])
        self.model = None
        
        #other parameters
        self.GPU_id = GPU_id

        #joint loss weight of the transfer model
        self.gamma = gamma

    
    '''
    :Section Model training of MultiSource_TargetLabel_dnnTL
    @Input： 
            X        expression matrix [samples x genes], pandas.core.frame.DataFrame object 
            y        0/1 label vector, pandas.core.frame.DataFrame object, whose row names correspond to the samples in data    
            model    [dictionary] of ClassifyModel objects used in training mode
    @Output：
            class of FX_M01_CF_M04_MTLP_TL_ChengMJ_BaseEstimator, where the model attribute stores the trained FX_M01_CF_M04_MTLP_TL_ChengMJ model
    '''           
    def fit(
            #fixed interface parameters
            self,  #this parameter can call the attributes and methods of the current class
            X: pd.DataFrame = None,    #training data matrix [must be named X], matrix format must be [samples x genes]
            y: pd.DataFrame = None,    #labels of training samples [must be named y]
            
            #other custom parameters
            model: Optional[dict] = None,  #input initialized source-model ClassifyModel collection
            GPU_id: Optional[int] = None   #which GPU to use, GPU_id can be 0,1,2,3 [None means CPU]
            ):
        '''
        :Part0 Data preprocessing''' 
        #(00) Device setting; by default use the first gpu ["cuda:0"], if using the second one set it as "cuda:1"
        GPU_id = GPU_id if GPU_id is not None else self.GPU_id  #if GPU is not specified in method parameters, search according to the class attribute
        device = torch.device(f"cuda:{GPU_id}" if GPU_id is not None else "cpu")
        data = pd.DataFrame(X).T
        label = pd.DataFrame(y)
        
        #(01) Handle multiple sources: rank sources
        if self.method_sort_sources == 'All':
            self.n_source=len(model)
        if self.method_sort_sources is not None and len(model)>1 and self.n_source < len(model): #if the number of sources is greater than 1 and a source-selection method is given
            if self.method_sort_sources == 'DTE':
                s_data_list = {name: t_model.data for name, t_model in model.items()}
                s_lab_list = {name: t_model.sample_label.values.flatten() for name, t_model in model.items()}
                sorted_sources = sort_source_by_DTE(s_data_list, s_lab_list, data, GPU_id=None)     #GPU is not used here because GPU memory is insufficient for large matrices; even so, it is 5 times faster than the previous function
                ordered_keys = sorted_sources.keys()
                # rank models according to source order
                model = OrderedDict((key, model[key]) for key in ordered_keys)
                                 
        #extract the first n_source source datasets for subsequent training
        ns = self.n_source
        source_data_list = [model[i].data.T.dropna(axis=1) for i in list(model.keys())[:ns]]   #prevent NA values in the matrix; if NA appears, delete NA gene columns with .dropna(axis=1)
        source_labels_list = [model[i].sample_label for i in list(model.keys())[:ns]]
        self.source_ls = list(model.keys())[:ns]  # record source list
        
        #(02) Handle target data
        T_data = data.T
        T_label = label
        T_data = T_data.dropna(axis=1)        #prevent NA values     
         
        #(03) Take the intersection genes of source and target data
        com_gene = T_data.columns
        for tmp_data in source_data_list[0:]:
            com_gene = com_gene.intersection(tmp_data.columns)
         
        source_data_list = [data[com_gene] for data in source_data_list]
        T_data = T_data[com_gene]
        #record input dimension and feature genes
        input_dim = T_data.shape[1]
        features= com_gene
        
        #(04) Append target data to the source-data list as the last element
        source_data_list.append(T_data)
        source_labels_list.append(T_label)

        
        '''
        :Part1 Data loading preparation''' 
        from collections import Counter
        # 1.1) Resample samples in each source domain to address class imbalance, and load them into data loaders 
        S_data_list_loaders = []
        for j in range(0,len(source_data_list)): # iterate through source_data_list and source_labels_list
            XTrain = source_data_list[j]
            YTrain = source_labels_list[j]
            # compute sample weights for each class
            class_sample_count = np.array([
                Counter(YTrain.iloc[:, 0])[0] / len(YTrain.iloc[:, 0]), 
                Counter(YTrain.iloc[:, 0])[1] / len(YTrain.iloc[:, 0])
            ])
        
            weight = 1. / class_sample_count  # compute class weights
            samples_weight = np.array([weight[t] for t in YTrain.values])  # assign sample weights
            samples_weight = torch.from_numpy(samples_weight)
            samples_weight = samples_weight.reshape(-1)  # Flatten out the weights
            # create a weighted random sampler
            sampler = WeightedRandomSampler(samples_weight.type('torch.DoubleTensor'), len(samples_weight), replacement=True)
            # data loading
            CDataset = TensorDataset(torch.FloatTensor(XTrain.values), torch.FloatTensor(YTrain.values.astype(int)))
            CLoader = DataLoader(dataset=CDataset, batch_size=self.batch_size, shuffle=False, sampler=sampler)  # use weighted random sampler
            S_data_list_loaders.append(CLoader) # store CLoader into the list
            
        # 1.2) Load the target domain into a data loader
        PDataset = torch.utils.data.TensorDataset(torch.FloatTensor(T_data.values), torch.FloatTensor(T_label.values.astype(int)) )
        PLoader = torch.utils.data.DataLoader(dataset=PDataset, batch_size=self.batch_size, shuffle=True)   #this data does not use labels [it is only used for mixing source domains and then classifying source vs target samples], so downsampling is unnecessary
        
        '''
        :Part2 Initialize neural network'''
        # 2.1) Initialize neural network and parameters
        DP = FX001_MTLP001_MTLP002(input_dim=input_dim, hiddenFX_dim=self.hiddenFX_dim, hiddenCF_dim=self.hiddenCF_dim, output_dim=self.output_dim, dropout_rate=self.dropout_rate)
        DP.to(device)
        
        # Ensure all parameters require gradients
        for param in DP.parameters():
            param.requires_grad = True
        
        # 2.2) Initialize parameter optimizer
        optimizer = torch.optim.Adagrad(itertools.chain(DP.parameters()), lr=self.lr)
        C_loss = torch.nn.BCELoss() # define the loss function (C_loss), using binary cross-entropy loss here
        
        '''
        :Part3 Train the model epoch by epoch'''
        cycled_iterators = [cycle(dat) for dat in S_data_list_loaders]  #put each source into an iterator separately
        for it in range(self.epoch):
            epoch_loss = []         
            '''
            :Point Feed each batch of samples into the model'''
            for batch_PLoader, *batch_S_data_list in zip(PLoader, *cycled_iterators): #extract data from each source domain and the target domain
                #set the model to training mode 
                DP.train()
                
                '''
                #(1) Import and extract the current batch of target data'''
                DataT = batch_PLoader
                xt = DataT[0].to(device)
                yt = DataT[1].view(-1, 1).to(device)
                
                '''
                #(2) Compute the current batch of target and source data separately; each source contributes the label-classifier loss loss_S_list, and each source together with the target contributes the domain-classifier loss loss_D_list'''
                loss_S_list = []
                loss_D_list = []
                ''' Introduced as a hyperparameter rather than set manually
                gamma is 𝛾 in the loss formula loss = 1/𝛾 𝑙𝑜𝑔⁡∑_(𝑖∈[1:𝑘])𝑒𝑥𝑝⁡(𝛾(𝑙𝑜𝑠𝑠_S_𝑖 + 𝑙𝑜𝑠𝑠_D_𝑚𝑎𝑥 )); larger values focus more on high-loss sources, while smaller values consider all losses more evenly
                @Practice shows: gamma should not be too large; if it is too large (e.g., >15), the loss may easily become Inf; it is best kept in the range 1-5'''
                gamma = self.gamma #1+ns*0.15         #set as a formula related to the number of sources ns; when there are more sources, the model focuses more on optimizing higher-loss sources, while with fewer sources it considers all losses more evenly
                for i_ns in range(1, len(batch_S_data_list)):  #for each source domain
                    '''
                    #(2.1) Import the current batch of source data and complete one forward pass'''
                    DataS = batch_S_data_list[i_ns]
                    xs = DataS[0].to(device)
                    ys = DataS[1].view(-1, 1).to(device)
                    if xs.size()[0] < 5 or xt.size()[0] < 5:  # skip training on the next batch if the number of samples in xs or xt is below a threshold (here, 5)
                        continue
                    #forward propagation, classify sample labels, and compute loss
                    yhat_xs = DP(xs)
                    loss1 = C_loss(yhat_xs, ys)
                    loss_S_list.append(loss1)
                    
                    ''' 
                    (2.2) Merge source-domain samples and target-domain samples, let the domain classifier determine whether a sample belongs to source or target, and compute the loss'''
                    Labels = torch.ones(xs.size(0), 1)  #set source-sample labels to 1
                    Labelt = torch.zeros(xt.size(0), 1) #set target-sample labels to 0
                    Lst = torch.cat([Labels, Labelt], 0).to(device)  #merge source and target labels
                    Xst = torch.cat([xs, xt], 0).to(device)          #merge source and target data
                    # classify source vs target samples and compute classification loss
                    yhat_DM = DP.dg_forward(Xst)  ##predicted domain label from global discriminator
                    loss2 = C_loss(yhat_DM, Lst)
                    loss_D_list.append(loss2)
        
                
                '''
                #(3) Aggregate the current batch losses from the label classifier and the domain classifier'''
                if len(loss_D_list)==0 or len(loss_S_list)==0: #if the sample size in this batch is <5, loss_S_list or loss_D_list will be empty, so skip this batch
                    continue
                #calculation must use torch built-in methods, otherwise gradients will be lost
                loss_D_max = torch.max(torch.stack(loss_D_list))  
                exp_sum = torch.sum(torch.exp(gamma * (torch.stack(loss_S_list) + loss_D_max)))  
                loss = (1 / gamma) * torch.log(exp_sum) 

                '''
                #(4) Backpropagate the loss and update model weights
                optimizer.zero_grad() #clear optimizer gradients to avoid accumulation: in each training iteration, gradients of model parameters need to be recomputed
                loss.backward() #backpropagation computes gradients; by the chain rule, gradients are propagated from the loss to every model parameter
                optimizer.step() #optimizer updates parameters using gradient information from the loss function to reduce the loss value [DP.fx_model.parameters(), DP.cf_model.parameters(), DM.cf_model.parameters() will be updated]
                '''
                # --- check whether loss is numerically abnormal ---
                if torch.isnan(loss) or torch.isinf(loss):
                    print(f"⚠️ Loss exploded at iter {it} ({loss.item()}), skipping the batch.")
                    continue  # skip this batch
                else:
                    # (4) Backpropagate the loss and update model weights
                    optimizer.zero_grad()  # clear optimizer gradients to avoid gradient accumulation
                    loss.backward()        # backpropagation computes gradients
                    optimizer.step()       # update model parameters

                # store batch-wise loss and performance metrics
                epoch_loss.append(loss)
        
            '''
            @Important After each batch is completed, delete the current variables on the GPU or detach them from the GPU
             Otherwise, repeated iterations using .to(device) may easily trigger the error: RuntimeError: CUDA error: device-side assert triggered '''
            del xt, yt, xs, ys, Lst, Xst
            #display the mean loss across batches in the current epoch
            loss_mean = torch.mean(torch.FloatTensor(epoch_loss))
            print("\n\nEpoch: {}".format(it))
            print("################loss_mean={}###############".format(loss_mean))
            
        #store and return the trained model
        DP.features = features
        DP.hyperparameters = {'input_dim':input_dim,
                               'hiddenFX_dim':self.hiddenFX_dim,
                               'hiddenCF_dim':self.hiddenCF_dim,
                               'dropout_rate':self.dropout_rate,
                               'output_dim':self.output_dim,
                               'batch_size':self.batch_size,
                               'epoch':self.epoch,
                               'lr':self.lr }
        #DP.scalerTrain = None #[multi-source transfer does not apply scaling, because there are multiple training matrices and it is unclear which matrix's scale parameters should be retained for new data] and [the referenced multi-source method Adversarial Multiple Source Domain Adaptation also does not involve scaling], so scaling is not used here
        DP.to(torch.device('cpu'))  #detach the model from GPU; otherwise, upon loading it will automatically occupy the GPU used during training
        self.model = DP
        del DP
        
        return(self)

    
    '''
    :Section Model prediction of FX_M01_CF_M04_MTLP_ChengMJ [#regular prediction function; output includes sample names]
    @Input： 
            X          expression matrix [samples X genes], pandas.core.frame.DataFrame object     
    @Output：
            y          predicted probability of positive samples
    '''
    def Predict(
                #fixed interface parameters
                self,
                X: pd.DataFrame = None,         #training data matrix [must be named X], matrix format must be [samples x genes]
                #other custom parameters
                GPU_id: Optional[int] = None):  #which GPU to use, GPU_id can be 0,1,2,3 [None means CPU]
        
        if self.model==None:
            raise RuntimeError("当前class需要训练神经网络模型\n")
        GPU_id = GPU_id if GPU_id is not None else self.GPU_id  #if GPU is not specified in method parameters, search according to the class attribute
        device = torch.device(f"cuda:{GPU_id}" if GPU_id is not None else "cpu")
        
        data = pd.DataFrame(X).T
        
        #extract features corresponding to the matrix
        features = self.model.features
        testX = data.T[features]
        
        #standardize the matrix using the trained scaler
        # if self.model.scalerTrain != None:
        #     testXN = self.model.scalerTrain.transform(testX.values)
        # else:
        testXN = testX.values
        testXN = torch.FloatTensor(testXN)
        
        #move the model to the specified device and switch to evaluation mode
        model = self.model.to(device)
        model.eval()
        
        #(1) Prediction
        yhat_list = []                 #store prediction results from each batch
        with torch.no_grad():              #gradient computation is usually unnecessary during prediction
            for i in range(0, testXN.size(0), self.batch_size):      #process data batch by batch using the batch size
                batch_testXN = testXN[i:i + self.batch_size]
                batch_testXN = batch_testXN.to(device)                         #move data to GPU device
                yhat_tmp = model(batch_testXN)                 #predict data
                yhat_list.append(yhat_tmp.cpu())
                print(i)
        yhat_temp = torch.cat(yhat_list, dim=0)                     #merge results from all batches
        
        #(2) Convert prediction-result format into two columns: CellID and Prediction
        yhat = yhat_temp.numpy()
        #convert data into a dictionary
        data_dict = {                                       
            "Label": testX.index.tolist(),
            "Prediction": yhat.flatten()  }
        #convert dictionary into DataFrame
        yhat_df = pd.DataFrame(data_dict)
        yhat_df.set_index("Label", inplace=True)       #set Label as the row index
        #delete variables stored on the device
        del batch_testXN, model
        
        return yhat_df


    ''' only used for the internal BaseEstimator interface
    @attention: the output of this method must be in array format, cannot carry sample names, and is only used for the interface'''
    def predict_proba(self, X,
                      GPU_id: Optional[int] = None):        #which GPU to use, GPU_id can be 0,1,2,3 [None means CPU]
        
        yhat_df = self.Predict(X=X, GPU_id=GPU_id)
        # create probability array: first column is negative-class probability, second column is positive-class probability
        t_prob = yhat_df['Prediction'].values
        proba = np.zeros((yhat_df.shape[0], 2))
        proba[:, 1] = t_prob
        proba[:, 0] = 1 - proba[:, 1]
        return proba  #probability array; first column is negative-class probability, second column is positive-class probability [the output data must follow this format]
    
  