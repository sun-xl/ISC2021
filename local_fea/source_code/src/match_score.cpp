#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include "imgcompare.h"
#include "sift_feature_extraction.h"
#include "sift_prefilter.h"
#include "sift_feature_match.h"

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


std::vector<std::pair<std::string, std::string> > readtxtfile(
    std::string vidlistpath) {
  std::ifstream frd(vidlistpath.c_str());
  std::vector<std::pair<std::string, std::string> > imgpair;
  while (!frd.eof()) {
    std::string line;
    frd >> line;
    if (line.size() == 0) {
      continue;
    }
    int npos = line.find(',');
    std::string src = line.substr(0, npos);
    std::string ref = line.substr(npos + 1);
    imgpair.push_back(std::make_pair(src, ref));
  }
  return imgpair;
}

int ImageSiftExtractor(std::string imgid, std::string imgpath, SIFT_FEATURE* siftfeat) {
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
  imageFeatureExtractionOne(dstframe, siftfeat, -1, 1);
  // 释放
  av_frame_free(&dstframe);

  if (imgframe != NULL) {
    av_frame_free(&imgframe);
  }

  return 0;
}

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << " <imgpair> <querydir> <refdir> <outfile>"
              << std::endl;
    return -1;
  }

  std::string imglistpath = argv[1];
  std::string srcdir = argv[2];
  std::string refdir = argv[3];
  std::string outfile = argv[4];
  std::vector<std::pair<std::string, std::string> > ilist = readtxtfile(imglistpath);

  std::ofstream OutFile(outfile.c_str());

  for (int i = 0; i < ilist.size(); i++) {
    std::string src = ilist[i].first;
    std::string ref = ilist[i].second;

    std::string src_path = srcdir + "/" + src;
    std::string ref_path = refdir + "/" + ref;

    // sift extract
    SIFT_FEATURE src_feat;
    SIFT_FEATURE ref_feat;
    int state = ImageSiftExtractor(src, src_path, &src_feat);
    int state2 = ImageSiftExtractor(ref, ref_path, &ref_feat);

    if (state != 0 || state2 != 0) {
      std::cout << src << "," << ref << "extract src or ref sift feat failed." << std::endl;
      delete_sift_feature_f(&src_feat);
      delete_sift_feature_f(&ref_feat);
      continue;
    }
    else {
      int matchKeyP = imgCompareTwo_xl(&ref_feat, &src_feat, 0, 1);
      delete_sift_feature_f(&src_feat);
      delete_sift_feature_f(&ref_feat);
      OutFile << src << "," << ref << ",score:" << matchKeyP << std::endl;
    }
  }
  OutFile.close();
  return 0;
}
