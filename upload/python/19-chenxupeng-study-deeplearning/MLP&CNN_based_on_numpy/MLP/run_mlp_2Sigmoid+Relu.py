from network import Network
from utils import LOG_INFO
from layers import Relu, Sigmoid, Linear
from loss import EuclideanLoss
from solve_net import train_net, test_net
from load_data import load_mnist_2d
import numpy as np


train_data, test_data, train_label, test_label = load_mnist_2d('data')

# Your model defintion here
# You should explore different model architecture
model = Network()
model.add(Linear('fc1', 784, 300, 0.01))
model.add(Sigmoid('Sigmoid1'))
model.add(Linear('fc2', 300, 300, 0.01))
model.add(Relu('Relu1'))
model.add(Linear('fc3', 300, 10, 0.01))

loss = EuclideanLoss(name='loss')

# Training configuration
# You should adjust these hyperparameters
# NOTE: one iteration means model forward-backwards one batch of samples.
#       one epoch means model has gone through all the training samples.
#       'disp_freq' denotes number of iterations in one epoch to display information.

config = {
    'learning_rate': 0.1,
    'weight_decay': 0.0001,
    'momentum': 0.7,
    'batch_size': 100,
    'max_epoch': 100,
    'disp_freq': 10,
    'test_epoch': 1
}
train_log_data = np.array(['']).astype('S')
test_log_data = np.array(['']).astype('S')
for epoch in range(config['max_epoch']):
    LOG_INFO('Training @ %d epoch...' % (epoch))
    train_net(model, loss, config, train_data, train_label, config['batch_size'], config['disp_freq'])
    log = train_net(model, loss, config, train_data, train_label, config['batch_size'], config['disp_freq'])
    train_log_data = np.concatenate((train_log_data,log))
    if epoch % config['test_epoch'] == 0:
        LOG_INFO('Testing @ %d epoch...' % (epoch))
        test_net(model, loss, test_data, test_label, config['batch_size'])
        log_ = test_net(model, loss, test_data, test_label, config['batch_size'])
        test_log_data = np.concatenate((test_log_data,log_))

np.savetxt('train_log',train_log_data,fmt='%s')
np.savetxt('test_log',test_log_data,fmt='%s')

