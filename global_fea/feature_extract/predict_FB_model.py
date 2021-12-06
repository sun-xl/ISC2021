import argparse
import os

import numpy as np
import torch
import h5py

from datasets.FB_dataset import TestFilelist
from torchvision import transforms as pth_transforms

from config import config
from config import update_config
from models.swin_transformer import get_cls_model, LinearClassifier, load_pretrained_weights


def build_argparse():
    parser = argparse.ArgumentParser(description='PyTorch CNN Fingerprint Net feature extract')
    parser.add_argument('--query', '-q', metavar='DATASET', default='datasets/query_list_crop.txt',
                        help='query dataset:  | (default: queryImgList.txt)')
    parser.add_argument('--total', '-t', metavar='DATASET', default='datasets/ref_list.txt',
                        help='total dataset:  | (default: totalImagList.txt)')
    parser.add_argument('--save_path', '-p',  default='h5_descriptors',
                        help="feature save path  (default: 'h5_descriptors')")
    parser.add_argument('--save_h5_name',  default='fb_descriptors.h5',
                        help="feature h5 file save path  (default: 'fb_descriptors.h5')")
    parser.add_argument('--model_type', '-mt',
                        default='EsViT_SwinB_W14',
                        help="assign your test model  (default: 'EsViT_SwinB_W14')")
    parser.add_argument('--model', '-m', default='checkpoints/final.pth',
                        help="the trained model path  ")
    parser.add_argument('--gpu-id', '-g', default='0', metavar='N',
                        help="gpu id used for testing (default: '0')")
    parser.add_argument('--workers', '-j', default=8, type=int, metavar='N',
                        help='number of data loading workers (default: 8)')
    parser.add_argument('--batch-size', '-b', default=250, type=int, metavar='N',
                        help='number of batch (default: 250)')
    args = parser.parse_args()
    return args


def main(args):
    # setting up the visible GPU
    os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu_id

    # loading network from path
    if args.model is None:
        print('please set the model file')

    print('using Model ' + args.model_type)
    print(">> Loading network:\n>>>> '{}'".format(args.model))

    # Init model
    if args.model_type == 'EsViT_SwinB_W14':
        # build model
        print('using EsViT_SwinB_W14')
        cfg_file = 'config/swin_base_patch4_window14_224.yaml'
        update_config(config, cfg_file)

        model = get_cls_model(config, is_teacher=True)

        swin_spec = config.MODEL.SPEC
        embed_dim=swin_spec['DIM_EMBED']
        depths=swin_spec['DEPTHS']
        n_last_blocks = 4
        model.n_last_blocks = n_last_blocks
        model.depths = depths
        num_features = []
        for i, d in enumerate(depths):
            num_features += [int(embed_dim * 2 ** i)] * d 
        num_features_linear = sum(num_features[-n_last_blocks:])

        model.head = LinearClassifier(num_features_linear, 256)

        # load parameters
        patch_size = 4
        load_pretrained_weights(model, args.model, "model", "swin_base", patch_size)
        net = model
    else:
        raise AttributeError('No %s model'%args.model_type)


    net.cuda()
    net.eval()

    transform = pth_transforms.Compose([
        pth_transforms.Resize((224, 224), interpolation=3),
        pth_transforms.ToTensor(),
        pth_transforms.Normalize((0.485, 0.456, 0.406), (0.229, 0.224, 0.225)),
    ])
    # creating dataset loader
    query_loader = torch.utils.data.DataLoader(
        TestFilelist(root='', flist=args.query, transform=transform),
        batch_size=args.batch_size, shuffle=False, num_workers=args.workers, pin_memory=True
    )
    query_data = extract_vectors(net, query_loader, "query_feature", args.save_path, save_txt=False)


    total_loader = torch.utils.data.DataLoader(
        TestFilelist(root='', flist=args.total, transform=transform),
        batch_size=args.batch_size, shuffle=False, num_workers=args.workers, pin_memory=True
    )
    ref_data = extract_vectors(net, total_loader, "features_feature", args.save_path, save_txt=False)

    
    # save in h5 file
    img_feat_h5 = os.path.join(args.save_path, args.save_h5_name)
    with h5py.File(img_feat_h5, "w") as ff:    #new H5 file
        qry_ids = ['Q' + str(x).zfill(5) for x in range(50_000)]
        ref_ids = ['R' + str(x).zfill(6) for x in range(1_000_000)]

        ff.create_dataset("query", data=query_data)
        ff.create_dataset("reference", data=ref_data)
        ff.create_dataset('query_ids', data=qry_ids)
        ff.create_dataset('reference_ids', data=ref_ids)


def extract_vectors(net, loader, version_name, save_path, save_txt=False):
    SaveTXT = True
    # move network to gpu and eval mode
    net.cuda()
    net.eval()
    #np.set_printoptions(precision=6)
    
    query_data_list = []

    ff = None
    if save_txt:
        save_dir = os.path.join(save_path, version_name)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        imgfeat = os.path.join(save_dir, 'imgfeat.txt')
        ff = open(imgfeat, 'w')

    all_batch = len(loader)
    # extracting vectors
    with torch.no_grad():
        for i, (input, name)in enumerate(loader):
            # for item in name: # save image list
            #     fn.write("%s\n" % item)
            input = input.cuda()
            out = net(input).cpu().data.squeeze().numpy()
            print('\r>>>> {}/{} done...'.format(i, all_batch), end='')

            # save in txt
            if save_txt:
                for i, item in enumerate(out):
                    img_name = name[i]
                    img_name = img_name[img_name.rfind('/')+1:img_name.rfind('.')]
                    ff.write(img_name + '||' + str(list(item))[1:-1]+'\n')
                    #ff.write("%s\n" % list(item))
                #print('\r>>>> {}/{} done...'.format(i, all_batch), end='')
            
            query_data_list.append(out)
    
    if save_txt:
        ff.close()
    
    data_np = np.vstack(query_data_list)

    print(data_np.shape)
    
    return data_np
    


if __name__=="__main__":
    args = build_argparse()
    main(args)
