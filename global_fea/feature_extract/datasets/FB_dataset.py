import os
import copy
import random

import torch
from PIL import Image

def default_loader(path):
    with open(path, 'rb') as f:
        return Image.open(f).convert('RGB')

def default_flist_reader(flist):
    """
    :param flist format: impath \nimpath \n ...
    :return: impath list
    """
    imlist = []
    with open(flist, 'r') as rf:
        for line in rf.readlines():
            impath = line.strip()
            imlist.append(impath)
    return imlist

class FBDataset(torch.utils.data.Dataset):
    def __init__(self, root, flist, transform_before = None, 
                 augmentation = None, transform_after = None,  
                 flist_reader=default_flist_reader, loader=default_loader):
        self.root = root
        self.imlist = flist_reader(flist)
        self.transform_before = transform_before
        self.augmentation = augmentation
        self.transform_after = transform_after
        self.loader = loader

    def __getitem__(self, index):
        img_name =  self.imlist[index]
        img = self.loader(os.path.join(self.root, img_name))
        if self.transform_before is not None:
            img = self.transform_before(img)
        # augmentation
        img_aug = copy.deepcopy(img)
        if self.augmentation is not None:
            img_aug = self.augmentation(img_aug)
        
        if self.transform_after is not None:
            img = self.transform_after(img)
            img_aug = self.transform_after(img_aug)

        input = torch.stack([img, img_aug], dim=0)
        target = torch.tensor([index, index])
        return input, target

    def __len__(self):
        return len(self.imlist)


def eval_datasets_flist_reader(flist):
    """
    :param flist format: impath,impath\nimpath impath\n ...
    :return: (impath, label) list
    """
    imlist = []
    with open(flist, 'r') as rf:
        for line in rf.readlines():
            q,r = line.strip().split(',')
            if q == 'query_id':
                continue
            imlist.append((q,r))
    return imlist



def eval_datasets_flist_reader_balance2(flist):
    """
    :param flist format: impath,impath\nimpath impath\n ...
    :return: (q, r, label) list
    """
    imlist = []
    label = 0
    with open(flist, 'r') as rf:
        for line in rf.readlines():
            q,r = line.strip().split(',')
            if q == 'query_id':
                continue
            imlist.append((q,r,label))
            label += 1
    # balance the labeled pairs and unlabeled pairs
    balance_list = []
    unlabeled_list = []
    for pair in imlist:
        q,r,_ = pair
        if r == '': # no ap pair
            unlabeled_list.append(pair)
        else:
            balance_list.append(pair)
    # 1:1
    balance_list.extend(random.sample(unlabeled_list, len(balance_list)))
    random.shuffle(balance_list)
    return balance_list


class FBEvalDataset_all_more_ref_adjust(torch.utils.data.Dataset):
    def __init__(self, root, flist, transform_before = None, 
                 augmentation = None, augmentation_query_single = None,
                 transform_after = None,  
                 flist_reader=eval_datasets_flist_reader_balance2, loader=default_loader,
                 ref_num=16):
        self.root = root
        self.imlist = flist_reader(flist)
        self.transform_before = transform_before
        self.augmentation = augmentation
        self.augmentation_query_single = augmentation_query_single
        self.transform_after = transform_after
        self.loader = loader
        self.ref_num = ref_num
        # 
        self.ref_labeled_map = dict()
        for i,(q,r,_) in enumerate(self.imlist):
            self.ref_labeled_map[r] = i         


    def __getitem__(self, index):
        q, r, label =  self.imlist[index]
        q = os.path.join(self.root, 'query_images', '%s.jpg'%q)

        if r == '': # no ap pair
            # r = "R%06d"%(random.randint(0,999999))
            img = self.loader(q)

            # augmentation
            img_aug = copy.deepcopy(img)
            if self.augmentation is not None:
                img_aug = self.augmentation(img_aug)
            if self.transform_after is not None:
                img = self.transform_after(img)
                img_aug = self.transform_after(img_aug)

            input = torch.stack([img_aug, img], dim=0)
            target = torch.tensor([label, label])
        else:
            r = os.path.join(self.root, 'reference_images', '%s.jpg'%r)
            
            img_q = self.loader(q)
            img_r = self.loader(r)

            if self.augmentation_query_single is not None:
                img_q = self.augmentation_query_single(img_q)

            if self.transform_after is not None:
                img_q = self.transform_after(img_q)
                img_r = self.transform_after(img_r)
            
            # Use (r, q) instead of (q, r), for using attacked image as 'anchor' in loss computation
            input = torch.stack([img_r, img_q], dim=0)  
            target = torch.tensor([label, label])
        
        # ref images, with no augmentation
        ref_num = self.ref_num
        img_nr_list = []
        ref_index_list = []
        for _ in range(ref_num):
            ref_index = random.randint(1,999999)
            neg_ref = "R%06d"%(ref_index)
            while (neg_ref in self.ref_labeled_map): # assign unlabeled ref
                ref_index = random.randint(1,999999)
                neg_ref = "R%06d"%(ref_index)
            img_nr = self.loader(os.path.join(self.root, 'reference_images', '%s.jpg'%neg_ref))
            if self.transform_after is not None: # normalization
                img_nr = self.transform_after(img_nr)

            img_nr_list.append(img_nr)
            ref_index_list.append([ref_index+30000])     
        img_nr_list = torch.stack(img_nr_list, dim=0)

        return input, target, img_nr_list, torch.tensor(ref_index_list)

    def __len__(self):
        return len(self.imlist)


class FBEvalDataset_all_more_ref(torch.utils.data.Dataset):
    def __init__(self, root, flist, transform_before = None, 
                 augmentation = None, augmentation_query_single = None,
                 transform_after = None,  
                 flist_reader=eval_datasets_flist_reader_balance2, loader=default_loader):
        self.root = root
        self.imlist = flist_reader(flist)
        self.transform_before = transform_before
        self.augmentation = augmentation
        self.augmentation_query_single = augmentation_query_single
        self.transform_after = transform_after
        self.loader = loader
        # 
        self.ref_labeled_map = dict()
        for i,(q,r,_) in enumerate(self.imlist):
            self.ref_labeled_map[r] = i         


    def __getitem__(self, index):
        q, r, label =  self.imlist[index]
        q = os.path.join(self.root, 'query_images', '%s.jpg'%q)

        if r == '': # no ap pair
            # r = "R%06d"%(random.randint(0,999999))
            img = self.loader(q)

            # augmentation
            img_aug = copy.deepcopy(img)
            if self.augmentation is not None:
                img_aug = self.augmentation(img_aug)
            if self.transform_after is not None:
                img = self.transform_after(img)
                img_aug = self.transform_after(img_aug)

            input = torch.stack([img_aug, img], dim=0)
            target = torch.tensor([label, label])
        else:
            r = os.path.join(self.root, 'reference_images', '%s.jpg'%r)
            
            img_q = self.loader(q)
            img_r = self.loader(r)

            if self.augmentation_query_single is not None:
                img_q = self.augmentation_query_single(img_q)

            if self.transform_after is not None:
                img_q = self.transform_after(img_q)
                img_r = self.transform_after(img_r)
            
            # Use (r, q) instead of (q, r), for using attacked image as 'anchor' in loss computation
            input = torch.stack([img_r, img_q], dim=0)  
            target = torch.tensor([label, label])
        
        # ref images, without augmentation
        ref_num = 16 #8
        img_nr_list = []
        ref_index_list = []
        for _ in range(ref_num):
            ref_index = random.randint(1,999999)
            neg_ref = "R%06d"%(ref_index)
            while (neg_ref in self.ref_labeled_map): # assign unlabeled ref
                ref_index = random.randint(1,999999)
                neg_ref = "R%06d"%(ref_index)
            img_nr = self.loader(os.path.join(self.root, 'reference_images', '%s.jpg'%neg_ref))
            if self.transform_after is not None: # normalization
                img_nr = self.transform_after(img_nr)

            img_nr_list.append(img_nr)
            ref_index_list.append([ref_index+30000])     
        img_nr_list = torch.stack(img_nr_list, dim=0)

        return input, target, img_nr_list, torch.tensor(ref_index_list)

    def __len__(self):
        return len(self.imlist)


class TestFilelist(torch.utils.data.Dataset):
    def __init__(self, root, flist, transform=None, 
                flist_reader=default_flist_reader, loader=default_loader):
        self.root = root
        self.imlist = flist_reader(flist)
        self.transfrom = transform
        self.loader = loader
        # print(len(self.imlist ))
        # print(self.imlist[:10])

    def __getitem__(self, index):
        """
        Args:
            index (int): Index
        Returns:
            tuple: (sample, target) where target is class_index of the target class.
        """
        img_name =  self.imlist[index]
        img = self.loader(os.path.join(self.root, img_name))

        if self.transfrom is not None:
            img = self.transfrom(img)

        return img, img_name

    def __len__(self):
        return len(self.imlist)

