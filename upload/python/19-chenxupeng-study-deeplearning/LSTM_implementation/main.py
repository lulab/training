#! /usr/bin/env python
import numpy as np
import tensorflow as tf
from tensorflow.python.framework import constant_op
import sys
import json
import time
import random
random.seed(1229)

from model import RNN, _START_VOCAB

tf.app.flags.DEFINE_boolean("is_train", True, "Set to False to inference.")
tf.app.flags.DEFINE_boolean("read_graph", False, "Set to False to build graph.")
tf.app.flags.DEFINE_integer("symbols", 18430, "vocabulary size.")
tf.app.flags.DEFINE_integer("labels", 5, "Number of labels.")
tf.app.flags.DEFINE_integer("epoch", 1e6, "Number of epoch.")
tf.app.flags.DEFINE_integer("embed_units", 300, "Size of word embedding.")
tf.app.flags.DEFINE_integer("units", 512, "Size of each model layer.")
tf.app.flags.DEFINE_integer("layers", 1, "Number of layers in the model.")
tf.app.flags.DEFINE_integer("batch_size", 16, "Batch size to use during training.")
tf.app.flags.DEFINE_string("data_dir", "../data", "Data directory")
tf.app.flags.DEFINE_string("train_dir", "./train", "Training directory.")
tf.app.flags.DEFINE_boolean("log_parameters", True, "Set to True to show the parameters")

FLAGS = tf.app.flags.FLAGS

def load_data(path, fname):
    print('Creating %s dataset...' % fname)
    data = []
    with open('%s/%s' % (path, fname)) as f:
        for idx, line in enumerate(f):
            tokens = line.split(' ')
            data.append({'label':tokens[0], 'text':tokens[1:]})
    return data

def build_vocab(path, data):
    print("Creating vocabulary...")
    vocab = {}
    for i, pair in enumerate(data):
        for token in pair['text']:
            if token in vocab:
                vocab[token] += 1
            else:
                vocab[token] = 1
    vocab_list = _START_VOCAB + sorted(vocab, key=vocab.get, reverse=True)
    if len(vocab_list) > FLAGS.symbols:
        vocab_list = vocab_list[:FLAGS.symbols]
    else:
        FLAGS.symbols = len(vocab_list)

    print("Loading word vectors...")
    #todo: load word vector from 'vector.txt' to embed, where the value of each line is the word vector of the word in vocab_list
    embed = np.zeros((len(vocab_list), 300), dtype='float32')
    word_index = {}
    for i, word in enumerate(vocab_list):
        word_index[word] = i
    with open('./data/vector.txt', 'r') as f:
        for line in f:
            c = line.strip().split()
            word = c[0]
            index = word_index.get(word)
            if index is not None:
                embed[index] = np.asarray(c[1:], dtype='float32')

    #embed = np.loadtxt('../data/vector.txt',dtype = 'str')[:,1:].astype(np.float32)
    return vocab_list, embed

def gen_batch_data(data):
    def padding(sent, l):
        return sent + ['_PAD'] * (l-len(sent))

    max_len = max([len(item['text']) for item in data])
    texts, texts_length, labels = [], [], []
        
    for item in data:
        texts.append(padding(item['text'], max_len))
        texts_length.append(len(item['text']))
        labels.append(int(item['label']))

    batched_data = {'texts': np.array(texts), 'texts_length':texts_length, 'labels':labels}

    return batched_data

def train(model, sess, dataset):
    st, ed, loss, accuracy = 0, 0, .0, .0
    gen_summary = True
    while ed < len(dataset):
        st, ed = ed, ed+FLAGS.batch_size if ed+FLAGS.batch_size < len(dataset) else len(dataset)
        batch_data = gen_batch_data(dataset[st:ed])
        outputs = model.train_step(sess, batch_data, summary=gen_summary)
        if gen_summary: 
            summary = outputs[-1]
            gen_summary = False
        loss += outputs[0]
        accuracy += outputs[1]
        sess.run(model.epoch_add_op)

    return loss / len(dataset), accuracy / len(dataset), summary

def evaluate(model, sess, dataset):
    st, ed, loss, accuracy = 0, 0, .0, .0
    gen_summary = True
    while ed < len(dataset):
        st, ed = ed, ed+FLAGS.batch_size if ed+FLAGS.batch_size < len(dataset) else len(dataset)
        batch_data = gen_batch_data(dataset[st:ed])
        outputs = sess.run([model.loss, model.accuracy, model.merged_summary_op],
                           {model.texts:batch_data['texts'],
                            model.texts_length:batch_data['texts_length'],
                            model.labels:batch_data['labels']})
        if gen_summary:
            summary = outputs[-1]
            gen_summary = False
        loss += outputs[0]
        accuracy += outputs[1]
    return loss / len(dataset), accuracy / len(dataset), summary

def inference(model, sess, dataset):
    st, ed, loss, accuracy = 0, 0, .0, .0
    gen_summary = True
    result = []
    while ed < len(dataset):
        st, ed = ed, ed+FLAGS.batch_size if ed+FLAGS.batch_size < len(dataset) else len(dataset)
        batch_data = gen_batch_data(dataset[st:ed])
        outputs = sess.run([model.predict_labels,
                            model.loss,
                            model.accuracy,
                            model.merged_summary_op],
                           {model.texts:batch_data['texts'],
                            model.texts_length:batch_data['texts_length'],
                            model.labels:batch_data['labels']})
        if gen_summary:
            summary = outputs[-1]
            gen_summary = False
        result += outputs[0].tolist()
        loss += outputs[1]
        accuracy += outputs[2]
    with open('result.txt', 'w') as f:
        for label in result:
            f.write('%d\n' % label)
    return loss / len(dataset), accuracy / len(dataset), summary

config = tf.ConfigProto()
config.gpu_options.allow_growth = True
with tf.Session(config=config) as sess:
    if FLAGS.is_train:
        print(FLAGS.__flags)
        data_train = load_data(FLAGS.data_dir, 'train.txt')
        data_dev = load_data(FLAGS.data_dir, 'dev.txt')
        data_test = load_data(FLAGS.data_dir, 'test.txt')
        vocab, embed = build_vocab(FLAGS.data_dir, data_train)
        
        model = RNN(
                FLAGS.symbols, 
                FLAGS.embed_units,
                FLAGS.units, 
                FLAGS.layers,
                FLAGS.labels,
                embed = embed,
                learning_rate=0.001)
        if FLAGS.log_parameters:
            model.print_parameters()
        
        if tf.train.get_checkpoint_state(FLAGS.train_dir):
            print("Reading model parameters from %s" % FLAGS.train_dir)
            model.saver.restore(sess, tf.train.latest_checkpoint(FLAGS.train_dir))
        else:
            print("Created model with fresh parameters.")
            tf.global_variables_initializer().run()
            op_in = model.symbol2index.insert(constant_op.constant(vocab),
                constant_op.constant(list(range(FLAGS.symbols)), dtype=tf.int64))
            sess.run(op_in)

        summary_writer = tf.summary.FileWriter('%s/log' % FLAGS.train_dir, sess.graph)
        while model.epoch.eval() < FLAGS.epoch:
            epoch = model.epoch.eval()
            random.shuffle(data_train)
            start_time = time.time()
            
            loss, accuracy, summary = train(model, sess, data_train)
            summary_writer.add_summary(summary, epoch)
            summary = tf.Summary()
            summary.value.add(tag='loss/train', simple_value=loss)
            summary.value.add(tag='accuracy/train', simple_value=accuracy)
            summary_writer.add_summary(summary, epoch)
            model.saver.save(sess, '%s/checkpoint' % FLAGS.train_dir, global_step=model.global_step)
            print("epoch %d learning rate %.4f epoch-time %.4f loss %.8f accuracy [%.8f]" % (epoch, model.learning_rate.eval(), time.time()-start_time, loss, accuracy))


            #todo: implement the tensorboard code recording the statistics of development and test set

            loss, accuracy, summary = evaluate(model, sess, data_dev)
            print("        dev_set, loss %.8f, accuracy [%.8f]" % (loss, accuracy))
            summary_writer.add_summary(summary, epoch)
            summary = tf.Summary()
            summary.value.add(tag='loss/dev', simple_value=loss)
            summary.value.add(tag='accuracy/dev', simple_value=accuracy)
            summary_writer.add_summary(summary, epoch)

            loss, accuracy, summary = inference(model, sess, data_test)
            print("        test_set, loss %.8f, accuracy [%.8f]" % (loss, accuracy))
            summary_writer.add_summary(summary, epoch)
            summary = tf.Summary()
            summary.value.add(tag='loss/test', simple_value=loss)
            summary.value.add(tag='accuracy/test', simple_value=accuracy)
            summary_writer.add_summary(summary, epoch)



