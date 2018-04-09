#! /usr/bin/env python
import numpy as np
import tensorflow as tf

from tensorflow.python.ops.nn import dynamic_rnn
from tensorflow.contrib.lookup.lookup_ops import MutableHashTable
from tensorflow.contrib.layers.python.layers import layers
from tensorflow.contrib.layers import l2_regularizer
from cell import GRUCell, BasicLSTMCell, BasicRNNCell
#from tensorflow.contrib.rnn import BasicRNNCell

PAD_ID = 0
UNK_ID = 1
_START_VOCAB = ['_PAD', '_UNK']

#W = tf.Variable(tf.truncated_normal([512, 5], stddev=0.1))
#B = tf.Variable(tf.constant(0.1, tf.float32, [5,]))

class RNN(object):
    def __init__(self,
            num_symbols,
            num_embed_units,
            num_units,
            num_layers,
            num_labels,
            embed,
            learning_rate=0.5,
            max_gradient_norm=5.0):
        #todo: implement placeholders
        self.texts = tf.placeholder(tf.string, [None, None], name='texts')  # shape: batch*len
        self.texts_length = tf.placeholder(tf.int32, [None], name='texts_length')  # shape: batch
        self.labels = tf.placeholder(tf.int32, [None], 'labels')  # shape: batch
        
        self.symbol2index = MutableHashTable(
                key_dtype=tf.string,
                value_dtype=tf.int64,
                default_value=UNK_ID,
                shared_name="in_table",
                name="in_table",
                checkpoint=True)
        # build the vocab table (string to index)
        # initialize the training process
        self.learning_rate = tf.Variable(float(learning_rate), trainable=False, dtype=tf.float32)
        self.global_step = tf.Variable(0, trainable=False)
        self.epoch = tf.Variable(0, trainable=False)
        self.epoch_add_op = self.epoch.assign(self.epoch + 1)


        self.index_input = self.symbol2index.lookup(self.texts)   # batch*len
        
        # build the embedding table (index to vector)
        if embed is None:
            # initialize the embedding randomly
            self.embed = tf.get_variable('embed', [num_symbols, num_embed_units], tf.float32)
        else:
            # initialize the embedding by pre-trained word vectors
            self.embed = tf.get_variable('embed', dtype=tf.float32, initializer=embed)
        
        
        self.embed_input = tf.nn.embedding_lookup(self.embed, self.index_input) #batch*len*embed_unit

        
        if num_layers == 1:
            cell = BasicLSTMCell(num_units)
        
        cell_dr = tf.nn.rnn_cell.DropoutWrapper(cell, input_keep_prob=0.5, output_keep_prob=0.5)
        outputs, states = dynamic_rnn(cell_dr, self.embed_input, self.texts_length, dtype=tf.float32, scope="rnn")
        labels = self.labels
        indices = tf.stack([tf.range(tf.shape(outputs)[0]),
                            self.texts_length - 1], axis=1)
        last_output = tf.gather_nd(outputs, indices)
        self.outputs = outputs
        self.states = states
        self.last_output = last_output
        #last_output = tf.reshape(last_output, (-1, num_units))
        logits = tf.layers.dense(last_output, 5)
        print(indices.shape)
        print(last_output.shape)



        self.logits = logits
        self.loss = tf.reduce_sum(tf.nn.sparse_softmax_cross_entropy_with_logits(labels=self.labels, logits=logits), name='loss')
        l2 = tf.cast(0.005,tf.float32) * sum(
                                             tf.nn.l2_loss(tf_var)
                                             for tf_var in tf.trainable_variables()
                                             if not ("noreg" in tf_var.name or "Bias" in tf_var.name)
                                             )
        self.loss += l2
        mean_loss = self.loss / tf.cast(tf.shape(self.labels)[0], dtype=tf.float32)
        predict_labels = tf.cast(tf.argmax(logits, 1, 'predict_labels'), tf.int32)
        self.predict_labels = predict_labels
        self.accuracy = tf.reduce_sum(tf.cast(tf.equal(self.labels, predict_labels), tf.int32), name='accuracy')

        self.params = tf.trainable_variables()
            
        # calculate the gradient of parameters
        #opt = tf.train.GradientDescentOptimizer(self.learning_rate)
        #opt = tf.train.RMSPropOptimizer(self.learning_rate)
        opt = tf.train.AdamOptimizer(self.learning_rate)
        gradients = tf.gradients(mean_loss, self.params)
        clipped_gradients, self.gradient_norm = tf.clip_by_global_norm(gradients, max_gradient_norm)
        self.update = opt.apply_gradients(zip(clipped_gradients, self.params), global_step=self.global_step)

        tf.summary.scalar('loss/step', self.loss)
        for each in tf.trainable_variables():
            tf.summary.histogram(each.name, each)
        tf.summary.histogram('logits', logits)
        tf.summary.histogram('gradient_norm', self.gradient_norm)
        tf.summary.histogram('rnn_output', self.last_output)
        for param in self.params:
            tf.summary.histogram('clipped_gradients/%s'%param.name, param)

        self.merged_summary_op = tf.summary.merge_all()
        
        self.saver = tf.train.Saver(write_version=tf.train.SaverDef.V2, 
                max_to_keep=3, pad_step_number=True, keep_checkpoint_every_n_hours=1.0)

    def print_parameters(self):
        for item in self.params:
            print('%s: %s' % (item.name, item.get_shape()))
    
    def train_step(self, session, data, summary=False):
        input_feed = {self.texts: data['texts'],
                self.texts_length: data['texts_length'],
                self.labels: data['labels']}
        output_feed = [self.loss, self.accuracy, self.gradient_norm, self.update]
        if summary:
            output_feed.append(self.merged_summary_op)
        return session.run(output_feed, input_feed)
