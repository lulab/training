from network import Network
from layers import Relu, Linear, Conv2D, AvgPool2D, Reshape
from utils import LOG_INFO
from loss import EuclideanLoss, SoftmaxCrossEntropyLoss
from solve_net import train_net, test_net
from load_data import load_mnist_4d
import numpy as np
import h5py

train_data, test_data, train_label, test_label = load_mnist_4d('../data')

# Your model defintion here
# You should explore different model architecture
model = Network()
model.add(Conv2D('conv1', 1, 4, 3, 1, 0.01))
model.add(Relu('relu1'))
model.add(AvgPool2D('pool1', 2, 0))  # output shape: N x 4 x 14 x 14
model.add(Conv2D('conv2', 4, 4, 3, 1, 0.01))
model.add(Relu('relu2'))
model.add(AvgPool2D('pool2', 2, 0))  # output shape: N x 4 x 7 x 7
model.add(Reshape('flatten', (-1, 196)))
model.add(Linear('fc3', 196, 10, 0.1))

loss = SoftmaxCrossEntropyLoss(name='loss')

# Training configuration
# You should adjust these hyperparameters
# NOTE: one iteration means model forward-backwards one batch of samples.
#       one epoch means model has gone through all the training samples.
#       'disp_freq' denotes number of iterations in one epoch to display information.

config = {
    'learning_rate': 0.12,
    'weight_decay': 0,
    'momentum': 0.78,
    'batch_size': 40,
    'max_epoch': 30,
    'disp_freq': 10,
    'test_epoch': 1
}


train_log_data = np.array(['']).astype('S')
test_log_data = np.array(['']).astype('S')
for epoch in range(config['max_epoch']):
    LOG_INFO('Training @ %d epoch...' % (epoch))
    train_net(model, loss, config, train_data, train_label, config['batch_size'], config['disp_freq'])
    log,output_1 = train_net(model, loss, config, train_data, train_label, config['batch_size'], config['disp_freq'])
    train_log_data = np.concatenate((train_log_data,log))
    if epoch % config['test_epoch'] == 0:
        LOG_INFO('Testing @ %d epoch...' % (epoch))
        test_net(model, loss, test_data, test_label, config['batch_size'])
        log_,output_2 = test_net(model, loss, test_data, test_label, config['batch_size'])
        test_log_data = np.concatenate((test_log_data,log_))

np.savetxt('train_log',train_log_data,fmt='%s')
np.savetxt('test_log',test_log_data,fmt='%s')
with h5py.File('image') as f:
    f.create_dataset('train_output',data = output_1)
    f.create_dataset('test_output',data = output_2)
