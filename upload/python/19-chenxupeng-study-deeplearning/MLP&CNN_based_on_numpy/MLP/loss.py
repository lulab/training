from __future__ import division
import numpy as np


class EuclideanLoss(object):
    def __init__(self, name):
        self.name = name

    def forward(self, input, target):
        '''Your codes here'''
        self.loss = 0.5* np.sum((input - target)**2)/input.shape[0]
        return self.loss

    def backward(self, input, target):
        '''Your codes here'''
        return (input - target)/input.shape[0]
