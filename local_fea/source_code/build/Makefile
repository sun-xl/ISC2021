CXX=g++
AR=ar
CP=cp
CXXPam=-c -g -std=c++0x

INCLUDE=-D_GLIBCXX_USE_CXX11_ABI=0 -I../extInc -I../extInc/opencv/ -I../util_inc/sift_inc -I../util_inc/fast_sift_inc  -Wno-sign-compare -Wno-write-strings  -Wno-unused-variable -Wno-strict-aliasing -Wno-return-type -Wno-deprecated-declarations -fpermissive -fopenmp
LDFLAGS=-Wl,-Bstatic ../extLib/libfast_sift_static.a -Wl,-Bdynamic ../extLib/libavcodec.so.57 ../extLib/libavdevice.so.57 ../extLib/libavfilter.so.6 ../extLib/libavformat.so.57 ../extLib/libavutil.so.55 ../extLib/libswresample.so.2 ../extLib/libswscale.so.4 ../extLib/libpostproc.so.54  ../extLib/libopencv_videoio.so.3.2 ../extLib/libopencv_imgproc.so.3.2 ../extLib/libopencv_core.so.3.2 ../extLib/libopencv_imgcodecs.so.3.2 ../extLib/libopencv_highgui.so.3.2 -Wl,-rpath=../extLib:../lib_cpu
SIFT_MODULE=../util_src/sift/

all:sift_cpu

sift_cpu:../src/siftfp_extract.cpp avf_conf.o util.o imgcompare.o featurematching.o
	@$(AR) rcs  ../lib_cpu/libsift_util.a avf_conf.o util.o imgcompare.o featurematching.o
	@$(CXX) -fPIC -O3 ../src/siftfp_extract.cpp -o localfea_extract_sift -Wl,-Bstatic ../lib_cpu/libsift_util.a $(LDFLAGS) $(INCLUDE) -Wl,-rpath=../lib_cpu:../extLib
	@$(CXX) -fPIC -O3 ../src/match_score.cpp -o match_score -Wl,-Bstatic ../lib_cpu/libsift_util.a $(LDFLAGS) $(INCLUDE) -Wl,-rpath=../lib_cpu:../extLib
	rm *.o

#SIFT OBJ
SIFT_OBJ: imgcompare.o featurematching.o
#	@echo -e $(COLORBG)"build sift obj"$(COLORFG)
featurematching.o: $(SIFT_MODULE)featurematching.cpp
	@$(CXX) -fPIC -c -std=c++0x $(INCLUDE) $(CXXPam) $(SIFT_MODULE)featurematching.cpp
imgcompare.o: $(SIFT_MODULE)imgcompare.cpp
	@$(CXX) -fPIC -c -std=c++0x $(INCLUDE) $(CXXPam) $(SIFT_MODULE)imgcompare.cpp

UTIL: avf_conf.o util.o
#	@echo -e "build cpu util"
avf_conf.o: ../util_src/sift/avf_conf.c ../util_inc/sift_inc/avf_conf.h
	@$(CC) -fPIC -c ../util_src/sift/avf_conf.c $(INCLUDE) $(LDFLAGS)
util.o: ../util_src/sift/util.c ../util_inc/sift_inc/util.h
	@$(CC) -fPIC -c ../util_src/sift/util.c $(INCLUDE) $(LDFLAGS)

clean:
	rm *.o
