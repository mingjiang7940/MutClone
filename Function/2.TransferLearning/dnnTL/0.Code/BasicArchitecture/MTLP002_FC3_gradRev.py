import torch
import torch.nn as nn
from torch.autograd import Function


'''
:Custom network propagation rule function: forward propagation remains unchanged, while backward propagation reverses the gradient
Implementation: during backpropagation, the model tries to minimize the difference between the positive and negative classes, rather than increasing the difference as in the conventional setting
'''
#####################################################
#Custom class GradReverse inherits from torch.autograd Function
class GradReverse(Function):
    # Static method: used for forward propagation
    @staticmethod
    def forward(ctx, x): #the main role of ctx is to pass information between forward propagation and backward propagation so that gradient computation can be performed in custom operations
        # 1. In forward propagation, directly return input x without any operation.
        return x.view_as(x)
        
    # Static method: used for backward propagation
    @staticmethod
    def backward(ctx, grad_output):
        # 2. In backward propagation, return a tensor with the same shape as grad_output,
        #    but with values equal to grad_output multiplied by -1.
        #    The purpose of this operation is to reverse the gradient direction.
        return (grad_output * -1)

#Custom function
def grad_reverse(x):
    return GradReverse.apply(x)  #For forward propagation of Function, apply should be called rather than forward
#######################################################




#Classify source samples (label = 1) and target samples (label = 0) based on shared features between source and target
# (Input) receives the extracted feature matrix [source samples and target samples x h_dim extracted features] for forward propagation
# Sequential neural network used to predict sample origin:
# (1) 1 Linear (fully connected) layer: performs weighted combination of features;
# (2) Dropout: during training, randomly sets some weights to 0 to prevent overfitting;
# (3) LeakyReLU: activation function that allows negative values to have small gradients
#     the above sequential unit is repeated 3 times
# (4) The final Linear layer has only one output value, used as the final prediction score
# (5) Sigmoid compresses the output value into the range 0-1, which is commonly used for classification problems and can be interpreted as the probability of belonging to the positive class
# | Linear | Dropout | LeakyReLU |
# | -----> | ------> | --------- |
# | Linear | Dropout | LeakyReLU |
# | -----> | ------> | --------- |
# | Linear | Dropout | LeakyReLU |
# | -----> | ------> | --------- |
# | Linear | Dropout | Sigmoid |
# | -----> | ------->|-------->|
# (Output) sample-origin score yhat
'''
:Neural network DG for predicting sample origin
@Input  merged matrix of source samples and target samples [rows are samples, columns are features extracted by FX]
@Output predicted score of sample origin (source samples are labeled 1, target samples are labeled 0)
'''
class MTLP002_FC3_gradRev(nn.Module):
    # Initialization method, accepts four parameters: input dimension (input_dim), hidden-layer dimension (hidden_dim), output dimension (hidden_dim), and dropout rate (dropout_rate)
    def __init__(self, input_dim, hidden_dim, output_dim=1, dropout_rate=0.5):
        super(MTLP002_FC3_gradRev, self).__init__()
        
        #Define a sequential fully connected neural network
        self.D1 = torch.nn.Sequential(
            nn.Linear(input_dim, hidden_dim),    # Linear layer from h_dim to h_dim
            nn.Dropout(p=dropout_rate), # Dropout layer
            nn.LeakyReLU(),             # Leaky ReLU activation function
            nn.Linear(hidden_dim, hidden_dim),    # Linear layer from h_dim to h_dim
            nn.Dropout(p=dropout_rate), # Dropout layer
            nn.LeakyReLU(),             # Leaky ReLU activation function
            nn.Linear(hidden_dim, hidden_dim),    # Linear layer from h_dim to h_dim
            nn.Dropout(p=dropout_rate), # Dropout layer
            nn.LeakyReLU(),             # Leaky ReLU activation function
            nn.Linear(hidden_dim, output_dim))        #output layer, reduces feature dimension to 1 as a single output value for sample prediction
        self.Drop1 = nn.Dropout(p=dropout_rate) #dropout layer

    #Forward propagation of the neural network using the custom method
    def forward(self, x):
        if x.shape[0]<=64*2:  #forward propagation is performed only for small batches [because the input here is the merged source-target matrix, the maximum batch size should be 2 times the maximum value, i.e., 2*64; otherwise the number of samples after merging in transfer learning will exceed the limit]
            x = grad_reverse(x)             #use custom function grad_reverse for forward propagation (standard forward propagation); gradients will be reversed when backward propagation is later called
            yhat = self.Drop1(self.D1(x))   #forward propagate to the output layer result and add a dropout layer
            return torch.sigmoid(yhat)      #use sigmoid function to compress the output value into the range 0-1
        else:
            print('batch size不要超过64\n')
            raise ValueError('Batch size too large. Maximum allowed is 64.')
        

