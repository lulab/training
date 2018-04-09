from __future__ import division
import numpy as np
from scipy import signal


class EuclideanLoss(object):
    def __init__(self, name):
        self.name = name

    def forward(self, input, target):
        return 0.5 * np.mean(np.sum(np.square(input - target), axis=1))

    def backward(self, input, target):
        return (input - target) / len(input)


class SoftmaxCrossEntropyLoss(object):
    def __init__(self, name):
        self.name = name

    def forward(self, input, target):

        a, c = input.shape
        exp = np.exp(input)
        ratio = exp / np.sum(exp, axis=1, keepdims=True)
        loss = -1.0 * np.sum(target * np.log(ratio)) / (a * c)
        self.ratio = ratio
        self.target, self.a, self.c =  target, a, c
        return loss

    def backward(self, input, target):

        target, a, c = self.target, self.a, self.c
        dX = self.ratio - target
        return dX / (a * c)
