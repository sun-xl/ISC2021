#include <iostream>
#include <fstream>
#include "sift_feature_extraction.h"
#include "imgcompare.h"
extern "C"
{
    #include "util.h"
}


int ComputeScaleRatio(int width, int height,
                                     float& scale_ratio) {
  scale_ratio = 1.0;
  int min_resize_pixel = 300;  // sunxl
  if (height > width) {
    scale_ratio = (min_resize_pixel * 1.0 / (width + 1e-5));
    if (scale_ratio >= 1.0) {
      scale_ratio = 1.0;
    }
  } else {
    scale_ratio = (min_resize_pixel * 1.0 / (height + 1e-5));
    if (scale_ratio >= 1.0) {
      scale_ratio = 1.0;
    }
  }
  return 0;
}


int main(int argc, char** argv) {
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " <imgid> <imgpath> <outfile>" << std::endl;
    return -1;
  }

  std::string imgid = argv[1];
  std::string imgpath = argv[2];
  std::string outfile = argv[3];
  std::ofstream OutFile(outfile.c_str());
  // 读取图片
  AVFrame* imgframe;
  int state = readImg((char*)imgpath.c_str(), &imgframe);
  if (state == 0) {
    std::cout << "read image error!" << std::endl;
    // av_frame_free(&imgframe);
    return -1;
  } 
  
  // resize
  float scale_ratio;
  if (imgframe->width <= 0 || imgframe->height <= 0) {
    return -1;
  }
  ComputeScaleRatio(imgframe->width, imgframe->height, scale_ratio);
  int resize_width = int(imgframe->width * scale_ratio);
  int resize_height = int(imgframe->height * scale_ratio);
  AVFrame* dstframe;
  initialFrame(&dstframe, imgframe->format, resize_width, resize_height);
  int scale_check;
  try {
    scale_check = scale_no_allocated(imgframe, dstframe, imgframe->format,
                                     resize_width, resize_height);
  }
  catch (std::exception& e) {
    av_frame_free(&dstframe);
    return -1;
  }
  
  // 提取
  SIFT_FEATURE siftfeat;
  imageFeatureExtractionOne(dstframe, &siftfeat, -1, 1);
 
  // 输出到文件
  int keypoint_num = siftfeat.keypoint_num;
  if (keypoint_num < 10) {
    return -1;
  }
  for (int i = 0; i < keypoint_num; ++i) {
    OutFile << imgid + "_" + std::to_string(i);
    for (int m = 0; m < 4; ++m) {
      OutFile << "_" << siftfeat.keypoint_ptr[4 * i + m];
    }
    OutFile << "||";
    for (int m = 0; m < 128; ++m) {
      OutFile << siftfeat.descriptor_ptr[128 * i + m];
      if (m < 128-1) {OutFile << ",";}
    }
    OutFile << std::endl;
  }
  OutFile.close();
  // 释放
  av_frame_free(&dstframe);
  delete_sift_feature_f(&siftfeat);

  if (imgframe != NULL) {
    av_frame_free(&imgframe);
  }

  return 0;
}
