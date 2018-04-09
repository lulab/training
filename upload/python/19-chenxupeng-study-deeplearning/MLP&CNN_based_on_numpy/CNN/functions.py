import numpy as np
from scipy import signal

def conv2d_forward(input, W, b, kernel_size, pad):

    n,c_in,h_in,w_in = input.shape
    padded_input = np.pad(input, pad_width=((0, 0), (0, 0), (pad, pad), (pad, pad)), mode='constant', constant_values=0)
    c_out = W.shape[0]
    h_out = h_in + 2*pad - kernel_size + 1
    w_out = w_in + 2*pad - kernel_size + 1
    output = np.zeros((n,c_out,h_out,w_out))
    for sample in range(n):
        for out_channel in range(c_out):
            for input_channel in range(c_in):
                output[sample,out_channel,:,:] +=signal.convolve2d(padded_input[sample,input_channel,:,:],np.rot90(W[out_channel,input_channel,:,:],2),mode='valid')
            output[sample,out_channel,:,:] += b[out_channel]
    return output



def conv2d_backward(input, grad_output, W, b, kernel_size, pad):

    n,c_in,h_in,w_in = input.shape
    n,c_out,h_out,w_out = grad_output.shape
    k = W.shape[3]
    grad_input = np.zeros((n,c_in,h_in+2*pad,w_in+2*pad))
    grad_W = np.zeros((c_out,c_in,k,k))
    grad_b = np.sum(np.sum(np.sum(grad_output,axis=3),axis=2),axis=0)
    for sample in range(n):
        for input_channel in range(c_in):
            for out_channel in range(c_out):
                grad_input[sample,input_channel,:,:]+=signal.convolve2d(grad_output[sample,out_channel,:,:],W[out_channel,input_channel,:,:],mode='full')
    grad_input = grad_input[:,:,pad:h_in+pad,pad:w_in+pad]
    for sample in range(n):
        for input_channel in range(c_in):
            for out_channel in range(c_out):
                grad_W[out_channel,input_channel,:,:] += signal.convolve2d(input[sample, input_channel, :, :], np.rot90(grad_output[sample, out_channel, :, :], 2), mode='valid')
    return grad_input,grad_W,grad_b


def avgpool2d_forward(input, kernel_size, pad):


    n, c_in, h_in, w_in = input.shape
    h_pad = h_in + 2 * pad
    w_pad = w_in + 2 * pad
    padded = np.pad(input, pad_width=((0, 0), (0, 0), (pad, pad), (pad, pad)), mode='constant', constant_values=0)

    h_out = int(h_pad / kernel_size)
    w_out = int(w_pad / kernel_size)

    output = np.zeros(shape=(n, c_in, h_out, w_out))
    col_num = h_pad * w_pad / (kernel_size * kernel_size)
    col = np.zeros(shape=(n, c_in, col_num, kernel_size * kernel_size))

    index = 0
    for i in range(h_out):
        for j in range(w_out):
            col[:, :, index, :] = \
            padded[:, :, i * kernel_size: (i + 1) * kernel_size, j * kernel_size: (j + 1) * kernel_size]\
                    .reshape(n, c_in, kernel_size * kernel_size)
            index += 1

    col_avg = np.average(col, axis=3)
    output = col_avg.reshape(n, c_in, h_out, w_out)
    return output

def avgpool2d_backward(input, grad_output, kernel_size, pad):

    n, c_in, h_in, w_in = input.shape
    _, _, h_out, w_out = grad_output.shape
    h_padded = h_in + 2 * pad
    w_padded = w_in + 2 * pad
    ones_matrix = np.ones((n, c_in, kernel_size, kernel_size))

    padded_grad_input = np.zeros((n, c_in, h_padded, w_padded))
    for i in range(h_out):
        for j in range(w_out):
            padded_grad_input[:, :, i * kernel_size: (i + 1) * kernel_size, j * kernel_size: (j + 1) * kernel_size] = np.multiply(ones_matrix, grad_output[:, :, i, j].reshape(n, c_in, 1, 1))
    grad_input =  padded_grad_input[:, :, pad: h_in + pad, pad: w_in + pad] / (kernel_size * kernel_size)
    return grad_input
