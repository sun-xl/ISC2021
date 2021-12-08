#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
#include <libavfilter/avfilter.h>

#ifdef __MINGW32__
#undef WIN32
#endif

#ifdef WIN32
_declspec(dllexport) void init();
_declspec(dllexport) int open_av(char filename[],AVFormatContext *fmt_ctx,AVCodecContext **vc_ctx,AVCodecContext **ac_ctx,AVCodec **vcodec,
								AVCodec **acodec,int *vi,int *ai,int *width,int *height,double *fps,double *duration);
_declspec(dllexport) int read_av(AVFormatContext *fmt_ctx,AVCodecContext *vc_ctx, AVCodecContext *ac_ctx, 
								AVCodec *vcodec, AVCodec *acodec,int vi, int ai, AVFrame **frame);
_declspec(dllexport) int set_avfilter(AVFilterGraph *filter_graph, AVFormatContext *fmt_ctx, AVStream *vstream, AVStream *astream, int sampleRate);
_declspec(dllexport) void print_stats_volumedetect(AVFilterContext *ctx);
_declspec(dllexport) int to_image(AVFrame* frame, char filename[]);
_declspec(dllexport) int to_image_png(AVFrame* frame, char filename[]);
_declspec(dllexport) int to_image_2(AVFrame* frame, char filename[], int width, int height, int crop, int encoder, int auto_gamma, double brightness, double contrast, double satuation);
_declspec(dllexport) int to_image_3(AVFrame* frame, char filename[], int width, int height, int crop, int encoder, int auto_gamma, double brightness, double contrast, double satuation, AVPacket* packet, int crf);
_declspec(dllexport) int yuv_image_adjust(AVFrame *frame, int auto_gamma, double brightness, double contrast, double satuation);
_declspec(dllexport) int hue_from_yuv(uint8_t y, uint8_t u, uint8_t v);
_declspec(dllexport) int detect_stripe(AVFrame* frame, int x, int y, int w, int h);
_declspec(dllexport) double in_the_ellipse(int y, int x, int cy, int cx, double sy, double sx, double theta);
_declspec(dllexport) void initialFrame(AVFrame** dest,int dest_fmt,int width, int height);
_declspec(dllexport) int scale_no_allocated(AVFrame* source,AVFrame* dest,int dest_fmt,int width, int height);
_declspec(dllexport) int scale(AVFrame* source,AVFrame** dest,int dest_fmt,int width, int height);
_declspec(dllexport) int readImg(char* filename,AVFrame** dest);
#endif

void init();
int open_av(char filename[],AVFormatContext *fmt_ctx,AVCodecContext **vc_ctx,AVCodecContext **ac_ctx,AVCodec **vcodec,
								AVCodec **acodec,int *vi,int *ai,int *width,int *height,double *fps,double *duration);
int open_av_headers(char filename[],AVFormatContext *fmt_ctx,AVCodecContext **vc_ctx,AVCodecContext **ac_ctx,AVCodec **vcodec,
								AVCodec **acodec,int *vi,int *ai,int *width,int *height, double *fps, double *duration,char *httpheaders);

int read_av(AVFormatContext *fmt_ctx,AVCodecContext *vc_ctx, AVCodecContext *ac_ctx, 
								AVCodec *vcodec, AVCodec *acodec,int vi, int ai, AVFrame **frame);
int read_av_keyFrame(AVFormatContext *fmt_ctx, AVCodecContext *vc_ctx, AVCodecContext *ac_ctx, 
								AVCodec *vcodec, AVCodec *acodec, int vi, int ai, AVFrame **frame);
int set_avfilter(AVFilterGraph *filter_graph, AVFormatContext *fmt_ctx, AVStream *vstream, AVStream *astream, int sampleRate);
void print_stats_volumedetect(AVFilterContext *ctx);

unsigned char **allocate_2d_uchar(int y, int x);
float **allocate_2d_float(int y, int x);
double **allocate_2d_double(int y, int x);
#ifdef WIN32
int round(double val);
#endif
unsigned char*** allocate_3d_uchar(int t, int x, int y);
double sqr(double x);
int to_image(AVFrame* frame, char filename[]);
int to_image_png(AVFrame* frame, char filename[]);
int to_image_2(AVFrame* frame, char filename[], int width, int height, int crop, int encoder, int auto_gamma, double brightness, double contrast, double satuation);
int to_image_3(AVFrame* frame, char filename[], int width, int height, int crop, int encoder, int auto_gamma, double brightness, double contrast, double satuation, AVPacket* packet, int crf);
int yuv_image_adjust(AVFrame *frame, int auto_gamma, double brightness, double contrast, double satuation);
int hue_from_yuv(uint8_t y, uint8_t u, uint8_t v);
int detect_stripe(AVFrame* frame, int x, int y, int w, int h);
double in_the_ellipse(int y, int x, int cy, int cx, double sy, double sx, double theta);
void initialFrame(AVFrame** dest,int dest_fmt,int width, int height);
int scale_no_allocated(AVFrame* source,AVFrame* dest,int dest_fmt,int width, int height);
int scale(AVFrame* source,AVFrame** dest,int dest_fmt,int width, int height);
int readImg(char* filename,AVFrame** dest);
int readImgFromBuffer(char* buffer, int buffer_size,AVFrame** dest);
