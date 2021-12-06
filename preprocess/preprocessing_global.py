import sys
import os
import cv2
import torch
from PIL import Image
import matplotlib.pyplot as plt
from shutil import copyfile




# Model init
yolov5_path = './yolov5'
pth_path =  './yolov5/runs/train/all_finetune3/weights/best.pt'
Thresh = 0.55

# load finetune weight
model = torch.hub.load(yolov5_path, 'custom', path=pth_path, source='local')


def select_max(results):
    
    max_rect = None
    results = results.xyxyn[0].cpu().numpy()
    
    max_conf = -1
    
    for result in results: 
        x0, y0, x1, y1, conf, _ = result
        
        if (y1-y0)*(x1-x0) > 0.85:
            continue
        if conf <= Thresh:
            continue
        
        ratio = (y1-y0)/(x1-x0)
        if ratio> 4 or ratio<1/4:
            continue
    
        if conf > max_conf:
            max_conf = conf
            max_rect = ( x0, y0, x1, y1 )
        
    return max_rect, max_conf


def pred_rect(img_path):
    img1 = Image.open(img_path)  # PIL image
    results = model(img1, size=640)  # includes NMS
    
    #挑选框
    width, height = img1.size
#     height, width = img.shape[:2]
    result, max_conf = select_max(results)
    
    if result:
        x0, y0, x1, y1 = result
        x0_pix, y0_pix, x1_pix, y1_pix = int(x0*width), int(y0*height), int(x1*width), int(y1*height)
        return (x0_pix, y0_pix, x1_pix, y1_pix), max_conf
    else:
        return None, None
    


# Load data

# img_dir = "../data/FB_competition/datasets/query_images"
# save_dir = "../data/FB_competition/datasets/query_images_chartlet_crop4"

img_dir = sys.argv[1]
save_dir = sys.argv[2]

img_list = list(os.listdir(img_dir))
if len(img_list) == 0:
    print('please input correct image folder path')
    exit()


i = 0
with open('global_crop_list.txt', 'w') as f:
    
    cur_img_list = img_list
    for img_path in cur_img_list:
        img_name = img_path
        print('%d/%d'%(i,len(cur_img_list)), img_name)
        img_path = os.path.join(img_dir, img_path)
        save_img_path = os.path.join(save_dir, img_name)
        i+=1

        square, max_conf = pred_rect(img_path)
        # show
        if square:
            img = cv2.imread(img_path)
#             cv2.rectangle(img, (square[0], square[1]), (square[2], square[3]), (255, 0, 0), thickness=5)
#             plt.imshow(img)
#             plt.show()

            crop_img = img[square[1]:square[3], square[0]:square[2],]
            cv2.imwrite(save_img_path, crop_img, [int(cv2.IMWRITE_JPEG_QUALITY), 100])
            f.write('%s %f\n' % (img_name, max_conf))
            print('valid rect detect !')
            
        else:
            #保存原图
#             img = cv2.imread(img_path)

#             #cv2.imwrite(save_img_path, img)
#             cv2.imwrite(save_img_path, img, [int(cv2.IMWRITE_JPEG_QUALITY), 100]) 
            copyfile(img_path, save_img_path)
            
            


