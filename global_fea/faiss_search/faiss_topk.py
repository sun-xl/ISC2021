import os
import os.path as osp
import numpy as np
import shutil
import sys
import argparse
import time
import faiss
import h5py


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
    print(imgfea.shape)

    return imgname, imgfea

if __name__=="__main__":
    in_feature = sys.argv[1]
    out_file = sys.argv[2]
    
    with h5py.File(in_feature, "r") as h5_file:
        query_imgfea = h5_file["query"][:]
        db_imgfea = h5_file["reference"][:]
        # Coerce IDs to native Python unicode string no matter what type they were before
        query_imgname = np.array(h5_file["query_ids"][:], dtype=object).astype(str).tolist()
        db_imgname = np.array(h5_file["reference_ids"][:], dtype=object).astype(str).tolist()
    
    dim = db_imgfea.shape[1]
    index = faiss.IndexFlatL2(dim)
    index.add(db_imgfea)

    f = open(out_file, 'w')

    # top one search
    D, I = index.search(query_imgfea, 1)
    for row,LIndex in enumerate(I):
        query_vid = query_imgname[row]
        for col,item in enumerate(LIndex):
            ref_vid = db_imgname[item]
            f.write(query_vid+'.jpg,'+ref_vid+'.jpg'+'\n')
            f.flush()

    f.close()
