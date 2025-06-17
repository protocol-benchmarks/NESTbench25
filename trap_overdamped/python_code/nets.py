
import torch


class FeedforwardNet(torch.nn.Module):
    """
    A feedforward neural network with layer normalization and Mish activation.

    This network maps input times to control protocol values.

    Attributes:
        model (torch.nn.Sequential): The sequential neural network model
    """

    def __init__(self, n_in, hidden_dim, num_layers, n_out):
        """
        Initialize the feedforward neural network.

        Args:
            n_in (int): Number of input features
            hidden_dim (int): Dimension of hidden layers
            num_layers (int): Number of hidden layers
            n_out (int): Number of output features
        """
        super().__init__()

        # Use Mish activation for better performance
        act = torch.nn.Mish()
        layers = []

        # Input layer
        layers.extend([
            torch.nn.Linear(n_in, hidden_dim),
            torch.nn.LayerNorm(hidden_dim, elementwise_affine=False),
            act
        ])

        # Hidden layers
        for _ in range(num_layers):
            layers.extend([
                torch.nn.Linear(hidden_dim, hidden_dim),
                torch.nn.LayerNorm(hidden_dim, elementwise_affine=False),
                act
            ])

        # Output layer (no activation or normalization)
        layers.append(torch.nn.Linear(hidden_dim, n_out))

        self.model = torch.nn.Sequential(*layers)

    def forward(self, x):
        """
        Forward pass through the network.

        Args:
            x (torch.Tensor): Input tensor

        Returns:
            torch.Tensor: Output tensor
        """
        return self.model(x)

    def save_model(self, file_path):
        """
        Save the model to a file.

        Args:
            file_path (str): Path to save the model
        """
        state_dict = {'model': self.model.state_dict()}
        torch.save(state_dict, file_path)

    def load_model(self, file_path):
        """
        Load the model from a file.

        Args:
            file_path (str): Path to the saved model
        """
        state_dict = torch.load(file_path)
        self.model.load_state_dict(state_dict['model'])
