# -*- coding: utf-8 -*-

import tensorflow as tf


def swish(x):
    return x * tf.nn.sigmoid(x)

def batchnorm(Ylogits, is_train, offset, convolutional=False):
    exp_moving_avg = tf.train.ExponentialMovingAverage(0.999) # adding the iteration prevents from averaging across non-existing iterations
    bnepsilon = 1e-5
    if convolutional:
        mean, variance = tf.nn.moments(Ylogits, [0, 1, 2])
    else:
        mean, variance = tf.nn.moments(Ylogits, [0])
    update_moving_averages = exp_moving_avg.apply([mean, variance])


    if is_train == True:
        m = exp_moving_avg.average(mean)
        v = exp_moving_avg.average(variance)

    else:
        m = mean
        v = variance
    Ybn = tf.nn.batch_normalization(Ylogits, m, v, offset, None, bnepsilon)
    return Ybn, update_moving_averages

def conv2d(x, filter_shape, strides_x, strides_y, padding, name):
    """
        Args:
        x: 4-D inputs. [batch_size, in_height, in_width, in_channels]
        filter_shape: A list of ints.[filter_height, filter_width, in_channels, out_channels]
        strides: A list of ints. 1-D tensor of length 4. The stride of the sliding window for each dimension of input.
        padding: A string from: "SAME", "VALID". The type of padding algorithm to use.
        Returns:
        h_conv:  A 4-D tensor. [batch_size, out_height, out_width, out_channels].
        if padding is 'SAME', then out_height==in_height.
        else, out_height = in_height - filter_height + 1.
        the same for out_width.
        """
    assert padding in ['SAME', 'VALID']
    strides=[1,strides_x, strides_y,1]
    with tf.variable_scope(name):
        W_conv = tf.get_variable('w', shape=filter_shape)
        b_conv = tf.get_variable('b', shape=[filter_shape[-1]])
        h_conv = tf.nn.conv2d(x, W_conv, strides=strides, padding=padding)
        h_conv_relu = tf.nn.relu(h_conv + b_conv)
    return h_conv_relu


def max_pooling(x, k_height, k_width, strides_x, strides_y, padding='SAME'):
    """max pooling layer."""
    ksize=[1,k_height, k_width,1]
    strides=[1,strides_x, strides_y,1]
    h_pool = tf.nn.max_pool(x, ksize, strides, padding)
    return h_pool

def dropout(x, keep_prob, name=None):
    """dropout layer"""
    return tf.nn.dropout(x, keep_prob, name=name)

def fc(x, in_size, out_size, name, activation=None):
    """fully-connect
        Args:
        x: 2-D tensor, [batch_size, in_size]
        in_size: the size of input tensor.
        out_size: the size of output tensor.
        activation: 'relu' or 'sigmoid' or 'tanh'.
        Returns:
        h_fc: 2-D tensor, [batch_size, out_size].
        """
    if activation is not None:
        assert activation in ['relu', 'sigmoid', 'tanh'], 'Wrong activation function.'
    with tf.variable_scope(name):
        w = tf.get_variable('w', shape = [in_size, out_size], dtype=tf.float32)
        b = tf.get_variable('b', [out_size], dtype=tf.float32)
        h_fc = tf.nn.xw_plus_b(x, w, b)
        if activation == 'relu':
            return tf.nn.relu(h_fc)
        elif activation == 'tanh':
            return tf.nn.tanh(h_fc)
        elif activation == 'sigmoid':
            return tf.nn.sigmoid(h_fc)
        else:
            return h_fc


K = 24  # first convolutional layer output depth
L = 48  # second convolutional layer output depth

W1 = tf.Variable(tf.truncated_normal([3, 3, 1, K], stddev=0.1))  # 6x6 patch, 1 input channel, K output channels
B1 = tf.Variable(tf.constant(0.1, tf.float32, [K]))
W2 = tf.Variable(tf.truncated_normal([3, 3, K, L], stddev=0.1))
B2 = tf.Variable(tf.constant(0.1, tf.float32, [L]))

W3 = tf.Variable(tf.truncated_normal([768, 10], stddev=0.1))
B3 = tf.Variable(tf.constant(0.1, tf.float32, [10]))




class Model:

    def __init__(self,
                 is_train,
                 learning_rate=0.001,
                 learning_rate_decay_factor=0.9995):
        self.x_ = tf.placeholder(tf.float32, [None, 1, 28, 28])
        self.y_ = tf.placeholder(tf.int32, [None])
        self.keep_prob = tf.placeholder(tf.float32)

        x = tf.reshape(self.x_, [-1, 28, 28, 1])

##############################################################################
##############################################################################

        # TODO: implement input -- Conv -- BN -- ReLU -- MaxPool -- Conv -- BN -- ReLU -- MaxPool -- Linear -- loss
        # The model
        # batch norm scaling is not useful with relus
        # batch norm offsets are used instead of biases
        stride = 1  # output is 28x28
        Y1l = tf.nn.conv2d(x, W1, strides=[1, stride, stride, 1], padding='SAME')
        Y1bn, update_ema1 = batchnorm(Y1l, is_train, B1, convolutional=True)
        Y1r = swish(Y1bn)
        h_pool1 = max_pooling(Y1r, 2, 2, 2, 2)
        Y1 = tf.nn.dropout(h_pool1,self.keep_prob)
        stride = 2  # output is 14x14
        Y2l = tf.nn.conv2d(Y1, W2, strides=[1, stride, stride, 1], padding='SAME')
        Y2bn, update_ema2 = batchnorm(Y2l, is_train,  B2, convolutional=True)
        Y2r = swish(Y2bn)
        h_pool2 = max_pooling(Y2r, 2, 2, 2, 2)
        Y2 = tf.nn.dropout(h_pool2,self.keep_prob)

        # reshape the output from the third convolution for the fully connected layer
        YY = tf.reshape(Y2, shape=[-1, 768])

        Ylogits = tf.matmul(YY, W3) + B3
        y_conv = tf.nn.softmax(Ylogits)

        self.update_ema = tf.group(update_ema1, update_ema2)

##############################################################################
##############################################################################

        self.loss = tf.reduce_mean(tf.nn.sparse_softmax_cross_entropy_with_logits(labels=self.y_, logits=y_conv))
        self.correct_pred = tf.equal(tf.cast(tf.argmax(y_conv, 1), tf.int32), self.y_)
        self.pred = tf.argmax(y_conv, 1)
        self.acc = tf.reduce_mean(tf.cast(self.correct_pred, tf.float32))

        self.learning_rate = tf.Variable(float(learning_rate), trainable=False,
                                         dtype=tf.float32)
        self.learning_rate_decay_op = self.learning_rate.assign(self.learning_rate * learning_rate_decay_factor)

        self.global_step = tf.Variable(0, trainable=False)
        self.params = tf.trainable_variables()
        self.train_op = tf.train.AdamOptimizer(self.learning_rate).minimize(self.loss, global_step=self.global_step,
                                                                            var_list=self.params)

        self.saver = tf.train.Saver(tf.global_variables(), write_version=tf.train.SaverDef.V2,
                                    max_to_keep=3, pad_step_number=True, keep_checkpoint_every_n_hours=1.0)


def weight_variable(shape):  # you can use this func to build new variables
    initial = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(initial)


def bias_variable(shape):  # you can use this func to build new variables
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial)

    # TODO: implemented the batch normalization func and applied it on conv and fully-connected layers
    # hint: you can add extra parameters (e.g., shape) if necessary

