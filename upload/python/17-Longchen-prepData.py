# -*- coding: utf-8 -*-
#!/usr/bin/env python
import random
import numpy as np
import matplotlib.pyplot as plt


def repeattimes(keylist):
    '''calculate occurring times of all tcrs'''
    times = {}
    for key in keylist:
        times[key] = times.get(key, 0) + 1
    return times


def keyfortime(dct, time):
    '''return keylist with occurring times of given num'''
    result = []
    keys = dct.keys()
    for key in keys:
        if dct[key] == time:
            result.append(key)
    return result


def keylengths(keylist):
    '''calculate lengths of tcrs in a set'''
    lengths = {}
    for i in keylist:
        key = len(i)
        lengths[key] = lengths.get(key, 0) + 1
    return lengths


def sample4freq(infile, time, key=0):
    '''return outfile with occurring time of time,
    according to keyth column array(keylist),
    e.g. 0-TCR, 5-antigen in all_merged.
    '''
    keylist = infile[:, key]
    dct = repeattimes(keylist)
    return_keylist = keyfortime(dct, time)
    mask = []
    # return_sample_num = return_keylist*time
    # return_samples = np.zeros((return_sample_num,infile.shape[1]))
    for k in keylist:
        if k in return_keylist:
            mask.append(True)
        else:
            mask.append(False)
    return_samples = infile[mask, :]
    return return_samples


def str2fix(key, length=18):
    '''turn TCR to fixed length, say 18, cut or add "#" which represent no amino axid;
    except for 0, which means return TCR without transforming'''
    if length == 0:
        return key
    else:
        n = len(key)
        if n >= length:
            new_key = key[:length]
        else:
            new_key = key + '#' * (length - n)
        return new_key


def aa2array(aminoAxid):
    '''change amino axid into a 20-dimension array, e.g. A-->(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)'''
    AAarray = np.zeros(20)
    AAloc = {}
    AAs = 'ARNDCQEGHILKMFPSTWYV'
    for i in range(len(AAs)):
        AA = AAs[i]
        AAloc[AA] = i
    if aminoAxid in AAs:
        loc = AAloc[aminoAxid]
        AAarray[loc] = 1
    return AAarray


def aas2mat(aminoAxidStr, length=18):
    '''change amino axid string (TCR) into a 20*18-dimension matrix'''
    new_aas = str2fix(aminoAxidStr, length=length)
    length = len(new_aas)
    AAsMat = np.zeros((20, length))
    for i in range(length):
        AA = new_aas[i]
        AAarray = aa2array(AA)
        AAsMat[:, i] = AAarray
    return AAsMat


def feature2mat(AAindexMap, aminoAxidStr, length=18, channel=0):
    '''add different AAindex feature array to given amino axid string (TCR), AAINDEX GRAPH NEEDED'''
    AAs = AAindexMap[0, 1:]                     # AAs = 'ARNDCQEGHILKMFPSTWYV'
    new_aas = str2fix(aminoAxidStr, length=length)
    length = len(new_aas)
    featureNum = AAindexMap.shape[0] - 2
    featureArray = np.zeros(length)
    featureMat = np.zeros((featureNum, length))
    featureValue = {}
    for i in range(featureNum):
        # with a 2-row header in AAindex Graph
        featureMap = AAindexMap[i + 2, 1:]
        for j in range(20):
            AA = AAs[j]
            featureValue[AA] = featureMap[j]
        for k in range(length):
            AA = new_aas[k]
            if AA == '#':
                featureArray[k] = 0
            else:
                featureArray[k] = featureValue[AA]
        featureMat[i] = featureArray
    featureMat = normalize(featureMat, channel=channel)
    return featureMat


def key2mat(AAindexMap, aminoAxidStr, length=18, channel=0):
    AAsMat = aas2mat(aminoAxidStr, length=length)
    featureMat = feature2mat(AAindexMap, aminoAxidStr,
                             length=length, channel=channel)
    keyMat = np.vstack((AAsMat, featureMat))
    return keyMat


def total2mat(AAindexMap, keys, head=0, length=18, mat_type='feature', channel=0):
    '''turn a list of samples into a matrix of feature, say f[:10] ==> 10*216
    mat_type: feature only, seq only, or both'''
    if length == 0:
        length = 18
    if mat_type == 'feature':
        sampleNum = len(keys)
        featureNum = AAindexMap.shape[0] - 2
        totalMat = np.zeros((sampleNum, featureNum * length))
        for i in range(sampleNum):
            featureMat = feature2mat(
                AAindexMap, keys[i][head:], length=length, channel=channel)
            totalMat[i, :] = featureMat.flatten()
    elif mat_type == 'seq':
        sampleNum = len(keys)
        totalMat = np.zeros((sampleNum, 20 * length))
        for i in range(sampleNum):
            aminoAxidStr = keys[i][head:]
            AAsMat = aas2mat(aminoAxidStr, length=length)
            totalMat[i, :] = AAsMat.flatten()
    elif mat_type == 'total':
        sampleNum = len(keys)
        featureNum = AAindexMap.shape[0] - 2
        totalMat = np.zeros((sampleNum, (featureNum + 20) * length))
        for i in range(sampleNum):
            aminoAxidStr = keys[i][head:]
            keyMat = key2mat(AAindexMap, aminoAxidStr,
                             length=length, channel=channel)
            totalMat[i, :] = keyMat.flatten()
    return totalMat


def sepData(AAindexMap, infile, head, TCR_col=0, Anti_col=6, part_Num_start=0, part_Num_end=10,
            ratio=0.7, length=18, mat_type='feature', return_type='part1', channel=0):
    keys = infile[:, TCR_col]
    labels = infile[:, Anti_col]
    sampleNum = len(keys)
    trainNum = int(ratio * sampleNum)
    train_mask = random.sample(range(sampleNum), trainNum)
    train_mask.sort()
    test_mask = np.delete(range(sampleNum), train_mask)
    train = infile[train_mask]
    test = infile[test_mask]
    np.set_printoptions(suppress=True)
    X_train = total2mat(
        AAindexMap, train[:, TCR_col], head=head, length=length, mat_type=mat_type, channel=channel)
    y_train = train[:, Anti_col]
    X_test = total2mat(
        AAindexMap, test[:, TCR_col], head=head, length=length, mat_type=mat_type, channel=channel)
    y_test = test[:, Anti_col]

    anti_times = repeattimes(labels)
    anti_freqs = list(set(anti_times.values()))
    anti_freqs.sort()
    anti_freqs.reverse()

    if return_type == 'part0':
        new_labels = np.zeros(labels.shape, dtype=int)
        for i in range(part_Num_start, part_Num_end):
            anti_new = keyfortime(anti_times, anti_freqs[i])
            new_labels += (labels == anti_new) * (i + 1)
        y_train = new_labels[train_mask]
        y_test = new_labels[test_mask]

    if return_type == 'part1':
        part_mask = np.zeros(labels.shape, dtype=int)
        for i in range(part_Num_start, part_Num_end):
            anti_new = keyfortime(anti_times, anti_freqs[i])
            part_mask += (labels == anti_new)
        part_mask = part_mask * np.arange(len(part_mask))
        # expel label 0; for ‘part2’, expel the 0th sample
        part_mask = list(set(part_mask))[1:]
        part_f1 = infile[part_mask]
        part_train_mask = random.sample(
            range(part_f1.shape[0]), int(ratio * len(part_mask)))
        part_train_mask.sort()
        part_test_mask = np.delete(range(part_f1.shape[0]), part_train_mask)

        part_train = part_f1[part_train_mask]
        part_test = part_f1[part_test_mask]
        np.set_printoptions(suppress=True)
        part_X_train = total2mat(
            AAindexMap, part_train[:, TCR_col], head=head, length=length, mat_type=mat_type, channel=channel)
        part_y_train = part_train[:, Anti_col]
        part_X_test = total2mat(
            AAindexMap, part_test[:, TCR_col], head=head, length=length, mat_type=mat_type, channel=channel)
        part_y_test = part_test[:, Anti_col]

        X_train = part_X_train
        X_test = part_X_test
        y_train = part_y_train
        y_test = part_y_test

    if return_type == 'part2':
        part_Num_start = part_Num_start
        part_Num_end = len(anti_freqs)
        X_train, y_train, X_test, y_test = \
            sepData(AAindexMap, infile, head, TCR_col=TCR_col, Anti_col=Anti_col,
                    part_Num_start=part_Num_start, part_Num_end=part_Num_end,
                    ratio=ratio, length=length, mat_type=mat_type, return_type='part1', channel=channel)

    return X_train, y_train, X_test, y_test


def replace(infile, Anti_col=6):
    anti_list = infile[:, Anti_col]
    anti_set = list(set(anti_list))
    replace_dict = {}
    keys = []
    k = 0
    for i in range(len(anti_set)):
        key = anti_set[i]
        if key not in keys:
            replace_dict[key] = k
            keys.append(key)
            k += 1
    return replace_dict


def anti2num(anti_list, infile, Anti_col=6):
    # anti_list = infile[:, Anti_col]
    num_list = []
    replace_dict = replace(infile, Anti_col=Anti_col)
    for anti in anti_list:
        num = replace_dict[anti]
        num_list.append(num)
    return np.asarray(num_list)


def normalize(mat, channel):  # channel = 0,1,2,3
    axis = 0
    if mat.min(axis=axis).all() == mat.max(axis=axis).all():
        pass
    else:
        channel_min_max = channel % 2
        channel_mean_std = channel // 2
        if channel_min_max == 1:
            mat = (mat - mat.min(axis=axis)) / \
                (mat.max(axis=axis) - mat.min(axis=axis))
        if channel_mean_std == 1:
            mat = (mat - mat.mean(axis=axis)) / (mat.std(axis=axis))
    return mat


def main():
    np.set_printoptions(suppress=True)
    f = np.loadtxt('TCR/Data/all_merged_GLIPH_input.txt', dtype='str')
    AAindexMap = np.loadtxt('TCR/Data/feature.txt', dtype=str)
    aminoAxidStr = 'VQGVPGPGVNEQFF'
    length = 0
    keys = f[:2, 0]

    AAsMat = aas2mat(aminoAxidStr, length=length)
    print(AAsMat)
    featureMat = feature2mat(AAindexMap, aminoAxidStr, length=length)
    print(featureMat)
    tcrMat = key2mat(AAindexMap, aminoAxidStr, length=length)
    print(tcrMat)
    
    totalMat = total2mat(AAindexMap, keys, head=1,
                         length=12, mat_type='feature')
    print(totalMat)
    print(totalMat.shape)
    print(np.sum(totalMat))


if __name__ == '__main__':
    main()
