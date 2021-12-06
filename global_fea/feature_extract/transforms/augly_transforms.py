import os
import random
import string

import augly.image as imaugs
import augly.utils as utils
from PIL import ImageFilter, Image
import torchvision
from torchvision import transforms as pth_transforms
import numpy as np

import augly.image.functional as F
import augly.image.utils as imutils
from typing import Any, Callable, Dict, List, Optional, Tuple, Union


random_RGB = lambda: (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))

string_pool = string.ascii_letters*2 + string.digits + ' '
string_pool_letter = string.ascii_letters + string.digits

def ramdom_string(length=10):
    letter_list = [random.choice(string_pool_letter)] + random.sample(string_pool, length-2) + [random.choice(string_pool_letter)]
    random_str = ''.join(letter_list)
    return random_str

class Meme(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        rgb = random_RGB()
        rgb_txt = random_RGB()
        text_len = random.randint(5, 9)
        text = ramdom_string(text_len)
        try:
            result_img = imaugs.meme_format(
                input_img,
                text=text,
                caption_height=random.randint(50, 200),
                meme_bg_color=rgb,
                text_color=rgb_txt,
            )
        except OSError:
            return input_img
        return result_img


class ShufflePixels(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        fact = random.randint(4, 6) * 0.1
        result_img = imaugs.shuffle_pixels(input_img, factor=fact)

        return result_img


class PixelizationRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        ratio = random.uniform(0.1, 0.5)
        result_img = imaugs.pixelization(input_img, ratio=ratio)

        return result_img


class BrightnessRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        factor = random.uniform(0.1, 2.0)
        result_img = imaugs.brightness(input_img, factor=factor)

        return result_img


class SaturationRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        factor = random.randint(1, 50) * 0.1
        result_img = imaugs.saturation(input_img, factor=factor)

        return result_img


class GrayscaleRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        gray_mode = random.choice(["luminosity", "average"])
        result_img = imaugs.grayscale(input_img, mode=gray_mode)

        return result_img


class BlurRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        factor = random.uniform(2, 10)
        result_img = imaugs.blur(input_img, radius=factor)

        return result_img


class SharpenRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        factor = random.randint(10, 30)
        result_img = imaugs.sharpen(input_img, factor=factor)

        return result_img

class JPEGEncodeAttackRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        q = random.randint(5, 20)
        result_img = imaugs.encoding_quality(input_img, quality=q)

        return result_img


class FilterRandom(object):
    def __init__(self):
        self.filter_list = [
            ImageFilter.MaxFilter,
            ImageFilter.UnsharpMask,
            #ImageFilter.CONTOUR,
            ImageFilter.EDGE_ENHANCE,
            ImageFilter.EDGE_ENHANCE_MORE,
            #ImageFilter.EMBOSS,
            ImageFilter.SMOOTH_MORE,
            ]

    def __call__(self, input_img):
        f = random.choice(self.filter_list)
        result_img = imaugs.apply_pil_filter(input_img, filter_type=f)

        return result_img


class PerspectiveTransform(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        sig = random.randint(10, 20) * 10.0
        aug = imaugs.PerspectiveTransform(sigma=sig)
        result_img = aug(input_img)

        return result_img


class OverlayOntoScreenshotScale(object):
    def __init__(self):
        self.screen_path = list(map(
            lambda x: os.path.join(utils.SCREENSHOT_TEMPLATES_DIR, x),
            ["web.png", "mobile.png", "mobile2.jpg", "web2.jpg", "mobile3.jpg",
            "mobile4.jpg", "mobile5.jpg", "mobile6.jpg", "web3.jpg"])
        )

    def __call__(self, input_img):
        aug = imaugs.Compose(
            [
                imaugs.OverlayOntoScreenshot(
                    template_filepath=random.choice(self.screen_path),
                    crop_src_to_fit=True,
                    resize_src_to_match_template = False,
                ),
                imaugs.Scale(factor=0.5),
            ]
        )
        result_img = aug(input_img)

        return result_img


class OverlayEmojiRandom(object):
    def __init__(self):
        self.emojis = os.path.dirname(utils.SMILEY_EMOJI_DIR,)
        self.emoji = []
        for root in os.listdir(self.emojis):
            root = os.path.join(self.emojis, root)
            emoji_list = []
            for file in os.listdir(root):
                if(file.lower().endswith(('.bmp', '.dib', '.png', '.jpg', '.jpeg', '.pbm', '.pgm', '.ppm', '.tif', '.tiff'))):
                    emoji_list.append(os.path.join(root, file))
            self.emoji.extend(emoji_list)

    def __call__(self, input_img):

        aug = imaugs.OverlayEmoji(
            emoji_path=random.choice(self.emoji),
            opacity=random.uniform(0.3, 1.0),
            emoji_size=random.uniform(0.2, 0.6),
            x_pos=random.randint(0, 80)*0.01,
            y_pos=random.randint(0, 80)*0.01,
        )
        result_img = aug(input_img)

        return result_img


def overlay_image2(
    image: Union[str, Image.Image],
    overlay: Union[str, Image.Image],
    output_path: Optional[str] = None,
    opacity: float = 1.0,
    overlay_size: float = 1.0,
    x_pos: float = 0.4,
    y_pos: float = 0.4,
    max_visible_opacity: float = 0.75,
    metadata: Optional[List[Dict[str, Any]]] = None,
    bboxes: Optional[List[Tuple]] = None,
    bbox_format: Optional[str] = None,
) -> Image.Image:
    assert 0.0 <= opacity <= 1.0, "Opacity must be a value in the range [0, 1]"
    assert 0.0 <= overlay_size <= 1.0, "Image size must be a value in the range [0, 1]"
    assert 0.0 <= x_pos <= 1.0, "x_pos must be a value in the range [0, 1]"
    assert 0.0 <= y_pos <= 1.0, "y_pos must be a value in the range [0, 1]"

    image = imutils.validate_and_load_image(image)

    func_kwargs = imutils.get_func_kwargs(metadata, locals())

    overlay = imutils.validate_and_load_image(overlay)

    im_width, im_height = image.size
    overlay_width, overlay_height = overlay.size
    new_height = max(1, int(im_height * overlay_size))
    new_width = int(overlay_width * new_height / overlay_height)
    overlay = overlay.resize((new_width, new_height))


    mask = Image.new(mode="L", size=overlay.size, color=int(opacity * 255))

    x = int(im_width * x_pos)
    y = int(im_height * y_pos)

    aug_image = image.convert(mode="RGBA")
    aug_image.paste(im=overlay, box=(x, y), mask=mask)

    imutils.get_metadata(
        metadata=metadata,
        function_name="overlay_image",
        aug_image=aug_image,
        **func_kwargs,
    )

    return imutils.ret_and_save_image(aug_image, output_path)


def overlay_emoji2(
    image: Union[str, Image.Image],
    output_path: Optional[str] = None,
    emoji_path: str = utils.EMOJI_PATH,
    opacity: float = 1.0,
    emoji_size: float = 0.15,
    x_pos: float = 0.4,
    y_pos: float = 0.8,
    metadata: Optional[List[Dict[str, Any]]] = None,
    bboxes: Optional[List[Tuple]] = None,
    bbox_format: Optional[str] = None,
) -> Image.Image:
    
    image = imutils.validate_and_load_image(image)

    func_kwargs = imutils.get_func_kwargs(metadata, locals())

    local_emoji_path = utils.pathmgr.get_local_path(emoji_path)

    aug_image = overlay_image2(
        image,
        overlay=local_emoji_path,
        output_path=output_path,
        opacity=opacity,
        overlay_size=emoji_size,
        x_pos=x_pos,
        y_pos=y_pos,
    )

    imutils.get_metadata(
        metadata=metadata,
        function_name="overlay_emoji",
        aug_image=aug_image,
        **func_kwargs,
    )

    return aug_image


class OverlayEmojiRandom2(object):
    def __init__(self):
        self.emojis = os.path.dirname(utils.SMILEY_EMOJI_DIR,)
        self.emoji = []
        for root in os.listdir(self.emojis):
            root = os.path.join(self.emojis, root)
            emoji_list = []
            for file in os.listdir(root):
                if(file.lower().endswith(('.bmp', '.dib', '.png', '.jpg', '.jpeg', '.pbm', '.pgm', '.ppm', '.tif', '.tiff'))):
                    emoji_list.append(os.path.join(root, file))
            self.emoji.extend(emoji_list)

    def __call__(self, input_img):
        
        result_img = overlay_emoji2(
                    input_img,
                    emoji_path=random.choice(self.emoji),
                    opacity=random.uniform(0.3, 1.0),
                    emoji_size=random.uniform(0.2, 0.7),
                    x_pos=random.randint(0, 50)*0.01,
                    y_pos=random.randint(0, 50)*0.01,
                )

        return result_img



class OverlayTextRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        
        text = []
        text_list = range(1000)
        width = random.randint(5,10)
        for _ in range(random.randint(1,3)):
            text.append(random.sample(text_list, width))
        
        text_size = random.uniform(0.1, 0.4)
            
        aug = imaugs.OverlayText(
            text=text,
            opacity=random.uniform(0.5, 1.0),
            font_size=random.uniform(0.1, 0.4),
            color=random_RGB(),
            x_pos=random.randint(0, 60) * 0.01,
            y_pos=random.randint(0, 60) * 0.01,
        )

        result_img = aug(input_img)

        return result_img


class PadRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        color = random_RGB()
        result_img = None
        if random.uniform(0, 1) > 0.5:
            result_img = imaugs.pad_square(
                input_img,
                color=color,
            )
        else:
            result_img = imaugs.pad(
                input_img,
                color=color,
                w_factor=random.uniform(0.1, 0.3),
                h_factor=random.uniform(0.1, 0.3),
            )

        return result_img


class OverlayRandom(object):
    def __init__(self):
        # '/group/20004/xinlongsun/joeqin/test_1130/test_pos/1245de22b70938bk.jpg'
        self.bg_list = []
        #with open('/group/20004/xinlongsun/joeqin/test_1130_list.txt', 'r') as rf:
        #with open('/cfs/cfs-c3ydde59/datasets/BG_img_list.txt', 'r') as rf:
        #    for line in rf.readlines():
        #        impath = line.strip()
        #        self.bg_list.append(impath)
        

    def __call__(self, input_img):

        bg = Image.open(random.choice(self.bg_list)).convert('RGB')

        width, height = bg.size  
        new_width = width * random.uniform(0.8, 1.0)
        if random.random() > 0.5:
            new_height = new_width * random.uniform(1.0, 1.5)
        else:
            new_height = new_width  * random.uniform(0.66, 1.0)
        left = (width - new_width)/2
        top = (height - new_height)/2
        right = (width + new_width)/2
        bottom = (height + new_height)/2
        bg = bg.crop((left, top, right, bottom))


        overlay_size = random.uniform(0.4, 0.8)

        max_len = max(input_img.size)

        y = random.uniform(0, max(1 - overlay_size - 0.1, 0.05))

        overlay_size_x = overlay_size*(new_width/new_height)*(input_img.size[0]/input_img.size[1])
        x = random.uniform(0, max(1 - overlay_size_x - 0.1, 0.05))

        result_img = imaugs.overlay_onto_background_image(input_img, 
                        background_image=bg,
                        x_pos=x, y_pos=y, opacity=1.0,
                        scale_bg=False, overlay_size=overlay_size)
        return result_img



#Add
class VerticalHorionalConvert(object):
    def __init__(self):
        self.wh_ratio = 9/16.
        #1280x720, 1138x640, 1024x576, 960x540,

    def __call__(self, img):
        w, h = img.size
        if w < h:   #vertical -> horional
            new_h = w * random.uniform(0.5, 1)
            h_ratio = new_h/h
            y1 = (1-h_ratio)/2
            y2 = y1 + h_ratio
            result_img = imaugs.crop(img, x1=0, x2=1, y1=y1, y2=y2)
                 
        else:   #horional -> vertical
            new_w = h * random.uniform(0.5, 1)
            w_ratio = new_w/w
            x1 = (1-w_ratio)/2
            x2 = x1 + w_ratio
            result_img = imaugs.crop(img, x1=x1, x2=x2, y1=0, y2=1)

        return result_img


class DouyinFilter(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        
        array_orig = np.array(input_img)
        
        if array_orig.shape[0] <= 20 or array_orig.shape[1] <=40:
            #print('[DouyinFilter]', array_orig.shape)
            return input_img
        
        array_r = np.copy(array_orig)
        array_r[:,:,1:3] = 255 #cyan
        
        array_b = np.copy(array_orig)
        array_b[:,:,0:2] = 255 #Y
        
        array_g = np.copy(array_orig)
        array_g[:,:,[0,2]] = 255 #R
 

        result_array = array_r[:-20, 40:, :] + array_b[20:, :-40, :] + array_g[10:-10, 20:-20, :]
        result_array[result_array>255] = 255
    
        result = Image.fromarray(result_array)

        return result


# https://arxiv.org/pdf/1805.09501.pdf, AutoAugment: Learning Augmentation Strategies from Data
class AutoAug(object):
    def __init__(self):
        self.policies = [pth_transforms.AutoAugmentPolicy.CIFAR10, pth_transforms.AutoAugmentPolicy.IMAGENET, pth_transforms.AutoAugmentPolicy.SVHN]

    def __call__(self, input_img):
        result = input_img
        for _ in range(2):  # multiple changes?
            policy = random.choice(self.policies)
            transform = pth_transforms.AutoAugment(policy) 
            result = transform(result)
            
        return result


class TorchvisionTrans(object):
    def __init__(self):
        self.trans_compose = [
            pth_transforms.RandomInvert(p=1), pth_transforms.RandomSolarize(threshold=192.0, p=1),
            pth_transforms.RandomEqualize(p=1), pth_transforms.RandomPosterize(bits=4, p=1), ]

    def __call__(self, input_img):
        transform = random.choice(self.trans_compose)
        result_img = transform(input_img)

        return result_img


class CropRandom(object):
    def __init__(self):
        pass

    def __call__(self, input_img):
        x1 = random.randint(10, 30) * 0.01
        y1 = random.randint(10, 30) * 0.01
        x2 = random.randint(70, 90) * 0.01
        y2 = random.randint(70, 90) * 0.01
        result_img = imaugs.crop(input_img, x1=x1, y1=y1, x2=x2, y2=y2)

        return result_img

color_jitter = torchvision.transforms.ColorJitter(brightness=0.4, contrast=0.4, saturation=0.2, hue=0.1)


augly_trans_list = [
    # overlay
    Meme(), OverlayOntoScreenshotScale(), OverlayEmojiRandom(), OverlayTextRandom(), OverlayEmojiRandom2(),
    # deformation
    PerspectiveTransform(), imaugs.RandomAspectRatio(), imaugs.RandomRotation(min_degrees=-180.0, max_degrees=180.0),
    PadRandom(), CropRandom(), OverlayRandom(), imaugs.HFlip(), VerticalHorionalConvert(), 
    pth_transforms.RandomAffine(degrees=(20, 80), translate=(0.2, 0.2), scale=(0.5, 0.9)),

    # noise
    ShufflePixels(), imaugs.random_noise, PixelizationRandom(), #imaugs.RandomPixelization(min_ratio = 0.05, max_ratio= 0.3), # 
    imaugs.Opacity(level=0.9), BrightnessRandom(), BlurRandom(),
    SharpenRandom(), JPEGEncodeAttackRandom(), FilterRandom(), color_jitter, SaturationRandom(),
    GrayscaleRandom(), DouyinFilter(), 
    pth_transforms.RandomInvert(p=1), pth_transforms.RandomPosterize(bits=4, p=1),
    pth_transforms.RandomSolarize(threshold=192.0, p=1), pth_transforms.RandomEqualize(p=1),  
    # Compose
    AutoAug(), # TorchvisionTrans(), 
]


if __name__ == '__main__':
    pass

