# -*- coding: utf-8 -*-
# @Time    : 2021-10-19 14:28
# @Author  : xinlongsun

import os
import os.path as osp
import numpy as np
import shutil
import sys
import argparse
import time
import faiss
import pickle

def load_feature(filepath):
    imgname = []
    imgfea = []
    with open(filepath, 'r') as f1:
        for line in f1.readlines():
            data = line.strip().split('|')
            imgname.append(data[0])
            img_featrue = np.array(eval(data[2]), dtype='float32')
            imgfea.append(img_featrue)
    imgfea = np.array(imgfea)

    return imgname, imgfea

if __name__=="__main__":
    in_dir = sys.argv[1]
    out_file = sys.argv[2]

    db_imgfea = pickle.load(open("./ref_sift_fea_300.pkl", "rb"))
    db_imgfea = db_imgfea.astype(np.float32)
    dim = db_imgfea.shape[1]

    cpu_index = faiss.IndexFlatL2(dim)
    co = faiss.GpuMultipleClonerOptions()
    co.shard = True
    co.useFloat16 = True
    index = faiss.index_cpu_to_all_gpus(cpu_index, co=co) # build the index
    index.add(db_imgfea)
    del db_imgfea
    
    db_imgname = pickle.load(open("./ref_sift_name_300.pkl", "rb"))

    f_out = open(out_file, 'w')

    for txt in os.listdir(in_dir):
        query_imgname, query_imgfea = load_feature(os.path.join(in_dir, txt))
        result_dict = dict()
        if len(query_imgname)==0:
            continue
        point_num = query_imgfea.shape[0]
        
        start = time.time()
        D, I = index.search(query_imgfea, 1)
        print(txt, ' faiss search time: ', time.time()-start)
        
        for row,LIndex in enumerate(I):
            query_vid = query_imgname[row]
            query_name = query_vid.split('_')[0]
            for col,item in enumerate(LIndex):
                ref_vid = db_imgname[item]
                ref_name = ref_vid.split('_')[0]
                if ref_name not in result_dict:
                    result_dict[ref_name] = [D[row,col]]
                else:
                    result_dict[ref_name].append(D[row,col])
        
        # sort result
        ref_list_new = sorted(result_dict.items(), key=lambda item: len(item[1]), reverse=True)

        if len(ref_list_new) > 1:
            if len(ref_list_new[0][1]) > 3 and len(ref_list_new[0][1]) > len(ref_list_new[1][1]):
                if np.array(ref_list_new[0][1]).mean() < 8000: # thre
                    f_out.write(query_name + '.jpg,' + ref_list_new[0][0] + '.jpg' + '\n')
                    f_out.flush()
        else:
            if len(ref_list_new[0][1]) > 5 and np.array(ref_list_new[0][1]).mean() < 8000:
                f_out.write(query_name + '.jpg,' + ref_list_new[0][0] + '.jpg' + '\n')
                f_out.flush()

    f_out.close()
