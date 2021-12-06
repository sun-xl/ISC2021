import random
import numbers
import os
from PIL import Image
from PIL import ImageFilter

import numpy as np

from . import functional as F

import cv2

import sys
import collections
if sys.version_info < (3, 3):
    Sequence = collections.Sequence
    Iterable = collections.Iterable
else:
    Sequence = collections.abc.Sequence
    Iterable = collections.abc.Iterable

root_path = os.getcwd()
#root_path = os.path.abspath(os.path.join(os.getcwd(), ".."))

__all__ = ["P_RandomCompose_FB", "P_Contrast", "P_hue", "P_Old", "P_Snow", "P_opencv_filter", "P_Noise", ]

transform_dict ={
    'P_RandomCompose_FB' : [1, 3],  
    'P_Contrast' : (0.5, 3),    
    'P_hue': (-0.5, 0.5),  
    'P_Old' : [1.0],  
    'P_Snow' : [1.0], 
    'P_Noise': [('salt', 0.01), ('gaussian', 10)], 
}


class P_RandomCompose_FB(object):
    """Random select n transforms and compose together.

	Args:
	    max transforms number: int

	Example:
	    >>>P_RandomCompose(transforms)
	"""

    def __init__(self, transforms, max_transform_num=None, min_transform_num=None):
        self.p = transform_dict['P_RandomCompose_FB'][0]
        self.transforms = transforms
        if max_transform_num == None:
            self.max_transform_num = transform_dict['P_RandomCompose_FB'][1]
        else:
            assert type(max_transform_num) == int
            self.max_transform_num = max_transform_num
        if min_transform_num == None:
            self.min_transform_num = transform_dict['P_RandomCompose_FB'][0]
        else:
            assert type(min_transform_num) == int
            self.min_transform_num = min_transform_num
        print("argumentation random min %d max %d :"%(self.min_transform_num, self.max_transform_num))

    def __call__(self, img):
        transform_num = random.randint(self.min_transform_num, self.max_transform_num)
        for t in random.sample(self.transforms, transform_num):
            try:
                img = t(img)
            except ValueError:
                print(t, 'error')
                img = img
            if img.mode == 'RGBA':
                img = img.convert('RGB')
        return img

    def __repr__(self):
        pass


class P_Contrast(object):
    """Randomly change the contrast of an image.

    Args:
        contrast (float or tuple of float (min, max)): How much to jitter contrast.
            contrast_factor is chosen uniformly from [max(0, 1 - contrast), 1 + contrast]
            or the given [min, max]. Should be non negative numbers.
    """
    def __init__(self):
        self.contrast = transform_dict['P_Contrast']
        self.contrast_factor = random.uniform(self.contrast[0], self.contrast[1])

    def __call__(self, img):
        return F.adjust_contrast(img, self.contrast_factor)

    def __repr__(self):
        format_string = self.__class__.__name__ + '('
        format_string += 'contrast={0})'.format(self.contrast_factor)
        return format_string


class P_Grayscale(object):
    def __init__(self):
        self.outchannel = transform_dict['P_Grayscale']

    def __call__(self, img):
        return F.to_grayscale(img, num_output_channels=self.outchannel)

    def __repr__(self):
        format_string = self.__class__.__name__ + '(outchannel={0})'.format(self.outchannel)
        return format_string


class P_opencv_filter(object):
    def __init__(self):
        self.colormaps = [cv2.COLORMAP_AUTUMN,
            cv2.COLORMAP_BONE,
            cv2.COLORMAP_COOL,
            cv2.COLORMAP_HOT,
            cv2.COLORMAP_HSV,
            cv2.COLORMAP_JET,
            cv2.COLORMAP_OCEAN,
            cv2.COLORMAP_PARULA,
            cv2.COLORMAP_PINK,
            cv2.COLORMAP_RAINBOW,
            cv2.COLORMAP_SPRING,
            cv2.COLORMAP_SUMMER,
            cv2.COLORMAP_WINTER]


    def __call__(self, img):
        #w, h = img.size

        cv_img = cv2.cvtColor(np.array(img), cv2.COLOR_RGB2GRAY)

        colormap = random.choice(self.colormaps)

        cv_img = cv2.applyColorMap(cv_img, colormap)

        return Image.fromarray(cv_img[:,:,::-1])

    def __repr__(self):
        return self.__class__.__name__ + '(opencv_filter)'


class P_hue(object):
    """Randomly change the hue of an image.

    Args:
        hue (float or tuple of float (min, max)): How much to jitter hue.
            hue_factor is chosen uniformly from [-hue, hue] or the given [min, max].
            Should have 0<= hue <= 0.5 or -0.5 <= min <= max <= 0.5.
    """

    def __init__(self):
        self.hue = transform_dict['P_hue']
        self.hue_factor = random.uniform(self.hue[0], self.hue[1])

    def __call__(self, img):
        return F.adjust_hue(img, self.hue_factor)

    def __repr__(self):
        format_string = self.__class__.__name__ + '('
        format_string += 'hue={0})'.format(self.hue_factor)
        return format_string


class P_Old(object):
    def __init__(self):
        weight = random.choice(transform_dict['P_Old'])
        if isinstance(weight, numbers.Number):
            self.weight = [float(weight), float(weight), float(weight)]
        else:
            self.weight = weight

    def __call__(self, img):
        if img.mode == 'L':
            print('Gray img cant add old filter, return origin img.')
            return img
        return F.filter_old(img, self.weight)

    def __repr__(self):
        return self.__class__.__name__ + '(weight={0})'.format(self.weight)


class P_Snow(object):
    def __init__(self):
        weight = random.choice(transform_dict['P_Snow'])
        assert isinstance(weight, float)
        self.weight = weight

    def __call__(self, img):
        if img.mode == 'L':
            print('Gray img cant add snow filter, return origin img.')
            return img
        return F.filter_snow(img, self.weight)

    def __repr__(self):
        return self.__class__.__name__ + '(weight={0})'.format(self.weight)


class P_Noise(object):
    def __init__(self):
        param = random.choice(transform_dict['P_Noise'])
        self.mode = param[0]
        self.ratio = param[1]

    def __call__(self, img):
        if img.mode == 'L':
            print('Gray img cant adFd noise filter, return origin img.')
            return img
        return F.filter_noise(img, self.mode, self.ratio)

    def __repr__(self):
        return self.__class__.__name__ + '(mode={0}, ratio={1})'.format(self.mode, self.ratio)

        