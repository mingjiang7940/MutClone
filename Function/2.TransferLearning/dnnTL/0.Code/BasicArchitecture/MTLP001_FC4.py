import torch.nn as nn
import torch.nn.functional as F

# Predict the response
# (1) Linear input layer: receives the Feature Extraction output, then applies BatchNorm1d, Dropout, and ReLU processing to prepare for downstream prediction
# | Linear(self.Sh) | BatchNorm1d(self.bn1) | Dropout(self.Drop) | ReLU
# | --------------> | --------------------> | -----------------> | ------> (ZX)
# (2) Sequential neural network for prediction:
# 2.1 Three Linear (fully connected) layers: perform weighted combinations of features;
# 2.2 Dropout: during training, randomly sets some weights to 0 to prevent overfitting;
# 2.3 LeakyReLU: activation function that allows negative values to have small gradients
# 2.4 The final Linear layer has only one output value, used as the final prediction score
# 2.5 Sigmoid compresses the output value into the range 0-1, which is commonly used for classification problems and can be interpreted as the probability of belonging to the positive class
# | Linear | Dropout | LeakyReLU |
# | ------------> | ------------> | -------------- |
# | Linear | Dropout | LeakyReLU |
# | ------------> | ------------> | -------------- |
# | Linear | Dropout | LeakyReLU |
# | ------------> | ------------> | -------------- |
# | Linear | Sigmoid |
# | ------------> | -------------- |
'''
:Note: the linear input layer F.relu(self.Drop(self.bn1(self.Sh((X))))) and
        nn.Sequential(nn.Linear(hidden_dim, hidden_dim),
                      nn.BatchNorm1d(hidden_dim),
                      nn.Dropout(p=dropout_rate),
                      nn.ReLU() )
        are two equivalent ways to define the network, and the latter is more convenient for fine-grained operations  '''
class MTLP001_FC4(nn.Module):
    # Initialization method, accepts four parameters: input dimension (input_dim), hidden-layer dimension (hidden_dim), output dimension (hidden_dim), and dropout rate (dropout_rate)
    ### drug response predictor ###
    def __init__(self, input_dim, hidden_dim, output_dim=1, dropout_rate=0.5):
        super(MTLP001_FC4, self).__init__()
        
        # Initialize layers and modules
        self.initial_layer = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),  # Linear layer from hidden_dim to hidden_dim
            nn.BatchNorm1d(hidden_dim),  # Batch normalization for hidden_dim
            nn.Dropout(p=dropout_rate)  # Dropout layer with specified dropout rate
        )
        
        self.Predictor = nn.Sequential(  # Sequential neural network
            nn.Linear(hidden_dim, hidden_dim),  # Linear layer from hidden_dim to hidden_dim
            nn.Dropout(p=dropout_rate),  # Dropout layer
            nn.LeakyReLU(),  # Leaky ReLU activation function
            nn.Linear(hidden_dim, hidden_dim),  # Linear layer
            nn.Dropout(p=dropout_rate),  # Dropout layer
            nn.LeakyReLU(),  # Leaky ReLU activation function
            nn.Linear(hidden_dim, hidden_dim),  # Linear layer
            nn.Dropout(p=dropout_rate),  # Dropout layer
            nn.LeakyReLU(),  # Leaky ReLU activation function
            nn.Linear(hidden_dim, output_dim),  # Linear layer from hidden_dim to 1
            nn.Sigmoid()  # Sigmoid activation function
        )
        
        # Initialize weights
        self._initialize_weights()

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='leaky_relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)

    ### Forward propagation, X is the input of the model
    def forward(self, X):
        ZX = self.initial_layer(X)  # Apply initial linear, batch norm, and dropout layers
        yhat = self.Predictor(ZX)  # Pass the result through the predictor network
        return yhat



        