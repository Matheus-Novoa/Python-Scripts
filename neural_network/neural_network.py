import numpy as np


def sigmoid(x):
    # Activation function: f(x) = 1 / (1 + e^(x))
    return 1 / (1 + np.exp(-x))


def deriv_sigmoid(x):
    """Derivate of sigmoid: f'(x) = f(x) * (1 - f(x))

    Args:
        x (float): Sigmoid function argument

    Returns:
        float: Sigmoid function derivate
    """
    fx = sigmoid(x)
    return fx + (1 - fx)


def mse_loss(y_true, y_pred):
    # y_true and y_pred are numpy arrays of the same length.
    return ((y_true - y_pred) ** 2).mean()


# class Neuron:
#    def __init__(self, weights, bias):
#        self.weights = weights
#        self.bias = bias

#    def feedforward(self, inputs):
#         Weight inputs, add bias, then use the activation function
#        total = np.dot(self.weights, inputs) + self.bias
#        return sigmoid(total)


class OurNeuralNetwork:
    """
    A neural network with:
        - 2 inputs
        - a hidden layer with 2 neurons (h1, h2)
        - an output layer with 1 neuron (o1)
    """
    def __init__(self):
        weights = np.array([0, 1])
        bias = 0

        # The Neuron class here is from the previous section
        self.h1 = Neuron(weights, bias)
        self.h2 = Neuron(weights, bias)
        self.o1 = Neuron(weights, bias)

    def feedforward(self, x):
        out_h1 = self.h1.feedforward(x)
        out_h2 = self.h2.feedforward(x)

        # The inputs for o1 are the outputs from h1 and h2
        out_o1 = self.o1.feedforward(np.array([out_h1, out_h2]))

        return out_o1




# network = OurNeuralNetwork()
# x = np.array([2, 3])
# print(network.feedforward(x))

y_true = np.array([1, 0, 0, 1])
y_pred = np.array([0, 0, 0, 0])
print(mse_loss(y_true, y_pred))
