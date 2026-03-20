import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Function
DevicePath = "D:/Project/0.MutClone/Function"
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/BasicArchitecture')
from FX001_FC1 import * 
from MTLP001_FC4 import * 
from MTLP002_FC3_gradRev import *
#torch.set_num_threads(1)  #number of CPU threads used by pytorch; this must be set, otherwise our server will default to using 52 threads

'''
:Assembled model      class of PredictFromExprs_TL for mutation prediction from expression using transfer learning
           
@Input  feature matrix of samples [samples x genes]
    
@Outpot prediction scores

@attention: [In practice, network training code is usually written as a separate function rather than placed inside the class defining the network]
#(1) The official PyTorch examples use separate training code rather than defining training functions inside the class. https://pytorch.org/tutorials/beginner/introyt/trainingyt.html
#(2) Placing training methods inside the class can easily cause many problems. https://discuss.pytorch.org/t/training-routine-inside-model-definition-class/24480 
#(3) Compared with placing training code into an independent training function, putting training code inside methods of the Network class definition has no clear advantage, and thousands of online examples have never adopted this design; therefore, it is best to use a separate training-function design.
#    https://jamesmccaffrey.wordpress.com/2023/06/07/why-isnt-pytorch-training-code-placed-inside-the-network-class-definition/
'''
class FX001_MTLP001_MTLP002(nn.Module):
   '''
   @attention Initialize the class with input_dim, hiddenFX_dim as the hidden dimension of the feature extractor, hiddenCF_dim as the hidden dimension of the classifier, output_dim=1, and dropout_rate=0.5'''
   def __init__(self, input_dim, hiddenFX_dim, hiddenCF_dim, output_dim=1, dropout_rate=0.5):
       super(FX001_MTLP001_MTLP002,self).__init__()
       
       #initialize submodels
       fx_model = FX001_FC1(input_dim=input_dim, hidden_dim=hiddenFX_dim, dropout_rate=dropout_rate) #FX_M01_CF_M04_MTLP_ChengMJ and FX_M01_CF_M03_MTLP_gradRev_ChengMJ share the same FX in forward propagation
       cf_model = MTLP001_FC4(input_dim=hiddenFX_dim, hidden_dim=hiddenCF_dim, output_dim=output_dim, dropout_rate=dropout_rate)
       dg_model = MTLP002_FC3_gradRev(input_dim=hiddenFX_dim, hidden_dim=hiddenFX_dim, output_dim=output_dim, dropout_rate=dropout_rate)
       
       #define three attributes of PredictFromExprs_TL for convenient extraction of submodels later
       self.fx_model=fx_model
       self.cf_model=cf_model
       self.dg_model=dg_model
       
   #  
   def forward(self, x):
       # default forward propagation
       F_x=self.fx_model(x)  #shared FX
       yhat_x=self.cf_model(F_x)
       return yhat_x
   
   #
   def dg_forward(self, x):
       #forward propagation when classifying samples as source domain or target domain; must be called explicitly
       F_x=self.fx_model(x)    #shared FX
       yhat_x=self.dg_model(F_x)
       return yhat_x
