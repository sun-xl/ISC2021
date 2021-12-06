from __future__ import division
import torch
import sys
import random
import math
from PIL import Image, ImageOps, ImageEnhance, ImageFilter, __version__
try:
    import accimage
except ImportError:
    accimage = None
import numpy as np
import numbers
import collections
import warnings

if sys.version_info < (3, 3):
    Sequence = collections.Sequence
    Iterable = collections.Iterable
else:
    Sequence = collections.abc.Sequence
    Iterable = collections.abc.Iterable


def adjust_contrast(img, contrast_factor):
    """Adjust contrast of an Image.

    Args:
        img (PIL Image): PIL Image to be adjusted.
        contrast_factor (float): How much to adjust the contrast. Can be any
            non negative number. 0 gives a solid gray image, 1 gives the
            original image while 2 increases the contrast by a factor of 2.

    Returns:
        PIL Image: Contrast adjusted image.
    """
    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    enhancer = ImageEnhance.Contrast(img)
    img = enhancer.enhance(contrast_factor)
    return img


def to_grayscale(img, num_output_channels=1):
    """Convert image to grayscale version of image.

    Args:
        img (PIL Image): Image to be converted to grayscale.

    Returns:
        PIL Image: Grayscale version of the image.
            if num_output_channels = 1 : returned image is single channel

            if num_output_channels = 3 : returned image is 3 channel with r = g = b
    """
    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    if num_output_channels == 1:
        img = img.convert('L')
    elif num_output_channels == 3:
        img = img.convert('L')
        np_img = np.array(img, dtype=np.uint8)
        np_img = np.dstack([np_img, np_img, np_img])
        img = Image.fromarray(np_img, 'RGB')
    else:
        raise ValueError('num_output_channels should be either 1 or 3')

    return img


def adjust_hue(img, hue_factor):
    """Adjust hue of an image.

    The image hue is adjusted by converting the image to HSV and
    cyclically shifting the intensities in the hue channel (H).
    The image is then converted back to original image mode.

    `hue_factor` is the amount of shift in H channel and must be in the
    interval `[-0.5, 0.5]`.

    See `Hue`_ for more details.

    .. _Hue: https://en.wikipedia.org/wiki/Hue

    Args:
        img (PIL Image): PIL Image to be adjusted.
        hue_factor (float):  How much to shift the hue channel. Should be in
            [-0.5, 0.5]. 0.5 and -0.5 give complete reversal of hue channel in
            HSV space in positive and negative direction respectively.
            0 means no shift. Therefore, both -0.5 and 0.5 will give an image
            with complementary colors while 0 gives the original image.

    Returns:
        PIL Image: Hue adjusted image.
    """
    if not(-0.5 <= hue_factor <= 0.5):
        raise ValueError('hue_factor is not in [-0.5, 0.5].'.format(hue_factor))

    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    input_mode = img.mode
    if input_mode in {'L', '1', 'I', 'F'}:
        return img

    h, s, v = img.convert('HSV').split()

    np_h = np.array(h, dtype=np.uint8)
    # uint8 addition take cares of rotation across boundaries
    with np.errstate(over='ignore'):
        np_h += np.uint8(hue_factor * 255)
    h = Image.fromarray(np_h, 'L')

    img = Image.merge('HSV', (h, s, v)).convert(input_mode)
    return img


def filter_old(img, channel_weight):
    """old filter to an Image.

    Args:
        img (PIL Image): PIL Image to be olded filter.
        channel_weight(float list): [r_weight, g_weight, b_weight] RGB channel's old weight.

    Returns:
        PIL Image: Old image.
    """
    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    img = np.array(img)
    rows, cols, dims = img.shape
    r_weight = channel_weight[0] + np.random.randint(-1000, 1000)/20000.0
    g_weight = channel_weight[1] + np.random.randint(-1000, 1000)/20000.0
    b_weight = channel_weight[2] + np.random.randint(-1000, 1000)/20000.0
    merge = np.empty([rows, cols, dims], float)
    merge[:, :, 0] = r_weight*(0.393* img[:, :, 0] + 0.769* img[:, :, 1] + 0.189* img[:, :, 2])
    merge[:, :, 1] = g_weight*(0.349* img[:, :, 0] + 0.686* img[:, :, 1] + 0.168* img[:, :, 2])
    merge[:, :, 2] = b_weight*(0.272* img[:, :, 0] + 0.534* img[:, :, 1] + 0.131* img[:, :, 2])

    img = np.uint8(np.clip(merge, 0, 255))
    img = Image.fromarray(img)

    return img


def filter_snow(img, weight):
    """Snow filter of an Image.

    Args:
        img (PIL Image): PIL Image channel to be snowed filter.
        weight(float): snow filter weight.

    Returns:
        PIL Image: Snow filter image.
    """
    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    maskR = np.array([0, 0, 0, 1, 1, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 10, 10, 11, 13, 13, 14, 14, 15, 16, 16, 17, 18, 19,
                  20, 21, 21, 22, 22, 23, 24, 26, 25, 27, 27, 29, 29, 30, 31, 33, 32, 34, 34, 36, 36, 37, 38, 39, 39,
                  40, 41, 43, 42, 43, 45, 47, 47, 48, 49, 51, 51, 52, 53, 55, 56, 57, 57, 58, 58, 58, 59, 61, 62, 63,
                  63, 65, 65, 66, 67, 69, 70, 71, 71, 73, 73, 74, 75, 77, 78, 79, 79, 81, 81, 82, 83, 85, 86, 87, 87,
                  89, 89, 90, 91, 93, 94, 95, 95, 97, 97, 98, 99, 101, 103, 104, 103, 106, 106, 107, 108, 110, 111, 112,
                  112, 114, 114, 115, 116, 118, 119, 119, 120, 122, 122, 123, 124, 126, 127, 128, 128, 129, 129, 130,
                  132, 134, 137, 138, 138, 140, 140, 139, 140, 142, 145, 145, 145, 147, 147, 149, 151, 153, 154, 154,
                  154, 156, 156, 157, 158, 160, 163, 164, 162, 164, 164, 167, 167, 170, 172, 172, 172, 174, 174, 176,
                  176, 178, 180, 180, 180, 182, 182, 183, 185, 187, 190, 191, 191, 191, 191, 192, 195, 197, 199, 200,
                  200, 202, 201, 202, 203, 207, 206, 207, 209, 209, 211, 213, 213, 216, 217, 217, 217, 220, 220, 220,
                  221, 224, 226, 227, 227, 230, 230, 231, 232, 234, 235, 235, 237, 238, 238, 240, 241, 243, 243, 246,
                  246, 248, 248, 249, 250, 253, 254, 255, 255], dtype=np.uint8)
    maskG = np.array([2, 2, 4, 7, 10, 12, 13, 13, 15, 15, 18, 18, 21, 22, 23, 24, 26, 26, 27, 29, 31, 32, 34, 33, 36, 36,
                  37, 39, 40, 41, 42, 42, 46, 46, 47, 48, 50, 51, 53, 53, 55, 55, 56, 57, 59, 60, 62, 62, 64, 64, 65,
                  66, 67, 69, 70, 69, 73, 72, 73, 75, 77, 79, 80, 79, 81, 81, 82, 83, 85, 86, 87, 87, 90, 90, 90, 91,
                  93, 94, 95, 95, 97, 97, 98, 99, 101, 102, 103, 103, 105, 105, 106, 107, 109, 110, 111, 111, 113, 113,
                  114, 115, 117, 118, 119, 119, 121, 121, 122, 123, 125, 126, 127, 127, 128, 128, 130, 131, 133, 133,
                  134, 135, 136, 136, 137, 138, 140, 141, 142, 142, 144, 144, 145, 146, 148, 149, 149, 148, 150, 150,
                  151, 152, 154, 155, 156, 156, 159, 159, 160, 160, 162, 163, 164, 164, 166, 166, 167, 168, 170, 171,
                  171, 171, 173, 173, 175, 175, 177, 178, 180, 180, 180, 180, 181, 182, 184, 184, 185, 186, 188, 188,
                  188, 188, 191, 192, 193, 193, 195, 195, 196, 196, 198, 198, 200, 200, 202, 202, 203, 203, 205, 206,
                  207, 207, 209, 209, 210, 211, 213, 213, 214, 214, 216, 217, 216, 217, 219, 220, 221, 221, 223, 223,
                  225, 225, 226, 227, 229, 229, 230, 230, 230, 231, 233, 235, 236, 236, 236, 236, 237, 238, 240, 241,
                  241, 241, 244, 242, 244, 245, 247, 247, 248, 248, 250, 250, 251, 252, 253, 254, 255, 255], dtype=np.uint8)
    maskB = np.array([1, 1, 3, 7, 9, 11, 12, 12, 16, 16, 18, 18, 21, 22, 24, 25, 26, 26, 27, 29, 31, 32, 33, 33, 35, 35, 38,
                  40, 41, 42, 43, 43, 46, 46, 47, 48, 50, 50, 52, 52, 54, 54, 55, 56, 60, 61, 63, 63, 65, 65, 66, 67,
                  68, 69, 70, 70, 73, 72, 73, 75, 77, 78, 79, 79, 81, 81, 82, 83, 85, 86, 87, 87, 89, 89, 89, 90, 92,
                  93, 94, 94, 96, 96, 97, 98, 100, 101, 102, 102, 104, 104, 105, 106, 108, 109, 110, 110, 112, 112, 113,
                  114, 116, 117, 118, 118, 120, 120, 121, 122, 124, 125, 126, 126, 130, 130, 129, 130, 132, 133, 134,
                  134, 136, 136, 137, 138, 140, 141, 142, 142, 144, 144, 145, 146, 148, 149, 149, 149, 151, 151, 152,
                  153, 155, 156, 157, 157, 159, 159, 160, 161, 163, 164, 165, 165, 167, 167, 168, 169, 171, 172, 172,
                  172, 172, 172, 174, 175, 177, 178, 179, 179, 180, 180, 181, 182, 184, 185, 186, 186, 188, 188, 189,
                  189, 192, 193, 194, 194, 196, 196, 195, 195, 197, 198, 199, 199, 201, 201, 202, 203, 205, 206, 207,
                  207, 209, 209, 210, 211, 213, 214, 215, 215, 217, 217, 216, 217, 219, 220, 221, 221, 223, 223, 225,
                  225, 227, 228, 229, 229, 231, 231, 231, 232, 232, 234, 235, 235, 236, 236, 237, 238, 240, 241, 241,
                  242, 244, 243, 245, 246, 248, 248, 247, 247, 249, 249, 250, 251, 253, 254, 255, 255], dtype=np.uint8)
    img = np.array(img)
    rows, cols, dims = img.shape
    snow_weight = weight + np.random.randint(-1000, 1000)/5000.0

    merge = np.empty([rows, cols, dims], np.int16)
    merge[:, :, 0] = maskR[img[:, :, 0]] * snow_weight
    merge[:, :, 1] = maskG[img[:, :, 1]] * snow_weight
    merge[:, :, 2] = maskB[img[:, :, 2]] * snow_weight

    img = np.uint8(np.clip(merge, 0, 255))
    img = Image.fromarray(img)
    return img


def filter_noise(img, noise_mode, noise_ratio=0.01):
    """Add noise to an Image.

    Args:
        img (PIL Image): PIL Image to be added noise.
        noise_mode (str): salt or gaussian noise.
        ratio (float): The ratio of noise strength. salt:0.01   gaussian:10

    Returns:
        PIL Image: Add noise image.
    """
    if not _is_pil_image(img):
        raise TypeError('img should be PIL Image. Got {}'.format(type(img)))

    assert noise_mode in ['salt', 'gaussian'], \
        'Padding mode should be either salt or gaussian'

    img = np.array(img)
    rows, cols, dims = img.shape
    if noise_mode == "salt":
        num = int(rows * cols * noise_ratio)
        for i in range(num):
            x = np.random.randint(0, rows)
            y = np.random.randint(0, cols)
            img[x, y, :] = 255

        for i in range(num):
            x = np.random.randint(0, rows)
            y = np.random.randint(0, cols)
            img[x, y, :] = 0

    if noise_mode == "gaussian":
        #img = img/255.0
        """
        img = skimage.util.random_noise(img, mode='gaussian', seed=None, clip=True)
        
        """
        mean = 0
        #var = 100
        sigma = noise_ratio
        gauss = np.random.normal(mean, sigma, (rows, cols, dims))
        gauss = gauss.reshape(rows, cols, dims)
        img = img + gauss

        img = np.uint8(np.clip(img, 0, 255))


    img = Image.fromarray(img)

    return img


