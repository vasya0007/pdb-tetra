import numpy as np
import json
import pprint


AMINO_ACIDS = ["A","R","N","D","B","C","E","Q","Z","G",
        "H","I","L","K","M","F","P","S","T","W","Y","V"]


def save_lables(data, file_name='lables.txt'):
    with open(file_name, 'w') as output:
        for id in data:
            output.write('\t'.join(id.split('_'))+'\n')


def cosine_distance(u, v):
    a = np.array(u)
    b = np.array(v)
    return 1 - np.dot(a,b)/np.sqrt(sum(a*a))/np.sqrt(sum(b*b))


def filter_duplicate(freq_dict):
    duplicate_ids = list()
    i=0
    for current_id in freq_dict:
        if current_id not in duplicate_ids:
            for another_id in freq_dict:
                if cosine_distance(freq_dict[current_id], freq_dict[another_id]) < 0.00001:
                    print(freq_dict[current_id])
                    print(freq_dict[another_id])
                    duplicate_ids.append(another_id)
        print('Номер Id:{} и число повторений:{}'.format(i, len(duplicate_ids)))
        i+=1
        #if i == 20:
        #    break
    print(len(freq_dict.keys()))
    for id in duplicate_ids:
        try:
            freq_dict.pop(id)
        except KeyError:
            pass
    print(len(freq_dict.keys()))
    return freq_dict


def create_freq_vector(seq):
    freq_vector = list()
    for amino_acid in AMINO_ACIDS:
        freq_vector.append(seq.count(amino_acid))
    return freq_vector


def filter_data(data):
    freq_dict = dict()
    for pdb_id in data:
        for chain in data[pdb_id]:
            if len(data[pdb_id][chain]) < 10:
                continue
            freq_dict[pdb_id+'_'+chain] = create_freq_vector(data[pdb_id][chain])
    return freq_dict


def open_data(path):
    with open(path, 'r') as raw_data:
        data = raw_data.read()
    return json.loads(data)


if __name__ == '__main__':
    chain_seq = open_data('../data_seq/chain_seq.json')
    
    freq_dict = filter_data(chain_seq)
    
    #print(len(freq_dict.keys()))
    treated_data = filter_duplicate(freq_dict)
    
    save_lables(treated_data)

    #print(filter_data(chain_seq))
    #pprint.pprint(filter_duplicate(freq_dict))
