import sys

result_dict = dict()

with open("global_pair_score.txt", "r") as f1:
    for line in f1.readlines():
        data = line.strip().split(',')
        img1 = data[0].split('.')[0]
        img2 = data[1].split('.')[0]
        sim = int(data[2].split(':')[-1])
        if img1+'_'+img2 not in result_dict:
            result_dict[img1+'_'+img2] = sim
        else:
            result_dict[img1+'_'+img2] += sim

with open("local_pair_score.txt", "r") as f1:
    for line in f1.readlines():
        data = line.strip().split(',')
        img1 = data[0].split('.')[0]
        img2 = data[1].split('.')[0]
        sim = int(data[2].split(':')[-1])
        if img1+'_'+img2 not in result_dict:
            result_dict[img1+'_'+img2] = sim
        else:
            result_dict[img1+'_'+img2] += sim

with open("crop_local_pair_score.txt", "r") as f1:
    for line in f1.readlines():
        data = line.strip().split(',')
        img1 = data[0].split('.')[0]
        img2 = data[1].split('.')[0]
        sim = int(data[2].split(':')[-1])
        if img1+'_'+img2 not in result_dict:
            result_dict[img1+'_'+img2] = sim
        else:
            result_dict[img1+'_'+img2] += sim

outfile = sys.argv[1]
f = open(outfile, 'w')
for key in result_dict:
    imgs = key.split('_')
    f.write(imgs[0]+'\t'+imgs[1]+'\t'+str(result_dict[key])+'\n')
    f.flush()
f.close()
