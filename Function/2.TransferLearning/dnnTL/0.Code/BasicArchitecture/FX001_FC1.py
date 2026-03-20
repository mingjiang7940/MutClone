import torch
import torch.nn as nn

# Feature extraction
# Feature Extraction: capture relevant information from the input feature vectors to achieve dimensionality reduction and denoising
# (input) accepts a [samp x genes] matrix for forward propagation
# The feature extractor architecture (self.EnE) includes:
# (1) 1 Linear (fully connected) layer: performs weighted combination of features;
# (2) BatchNorm layer: batch normalization applied to the Linear layer output;
# (3) ReLU: activation function that sets negative values to 0;
# (4) Dropout: during training, randomly sets some weights to 0 to prevent overfitting;
# The final output layer contains hidden_dim extracted features
# | Linear | BatchNorm1d | ReLU | Dropout
# | -----> | ----------> | ---> | ------> 
# (output) outputs a [samp x hidden_dim] matrix as the feature extraction result
'''
:Part1_FeatureExtractor_FX
@Input  feature matrix of samples [samples x features]
@Output extracted feature matrix [samples x hidden-layer features]
'''
class FX001_FC1(nn.Module):
    # Initialization method, accepts four parameters: input dimension (input_dim), hidden-layer dimension (hidden_dim), output dimension (output_dim), and dropout rate (dropout_rate)
    def __init__(self, input_dim, hidden_dim, dropout_rate=0.5):
        # Call the parent class nn.Module initialization method to ensure all parent attributes are properly initialized
        super(FX001_FC1, self).__init__()
        
        # Define a neural network module (Sequential) containing multiple layers
        self.EnE = nn.Sequential(
            nn.Linear(input_dim, hidden_dim), # Linear layer, transforms input dimension input_dim to hidden dimension hidden_dim
            nn.BatchNorm1d(hidden_dim),       # Batch normalization layer, normalizes the output of the hidden layer
            nn.ReLU(),                        # ReLU activation function to introduce non-linearity
            nn.Dropout(p=dropout_rate)        # Dropout layer, randomly drops some neurons during training to reduce overfitting
        )
    
    # Define the forward propagation method of the neural network, where x is the model input
    def forward(self, x):
        output = self.EnE(x) # Pass input x through the defined neural network module
        return output        # Return the output of the module
    