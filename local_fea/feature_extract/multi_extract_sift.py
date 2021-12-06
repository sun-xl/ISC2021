import os
from concurrent.futures import ProcessPoolExecutor
import time
import sys

root_path = sys.argv[2]
out_path = sys.argv[3]

def multi_download(rowkey_list):
    imgid = rowkey_list[0]
    imgpath = rowkey_list[1]
    outfile = os.path.join(out_path,imgid+'.txt')
    try:
        os.system('./localfea_extract_sift {} {} {}'.format(imgid, imgpath, outfile))
    except Exception as e:
        print('error:'+'\t'+imgid+'\t'+imgpath)


def main(image_list):
    start = time.time()
    download_list = []

    with open(image_list, 'r') as f1:
        for line in f1.readlines():
            imgname = line.strip().split(' ')[0]
            imgid = imgname.split('.')[0]
            imgpath = os.path.join(root_path, imgname)
            download_list.append([imgid, imgpath])
    print(len(download_list))
     
    with ProcessPoolExecutor(max_workers=100) as executor:
        try:
            executor.map(multi_download, download_list)
        except Exception as e:
            print(e)
    print(time.time()-start) 

if __name__ == "__main__":
    image_list = sys.argv[1]
    main(image_list)
    
