//#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
#include <libavfilter/avfilter.h>
#include "util.h"
#include "avf_conf.h"
//#include <iostream>
//#include <glog/logging.h>
//#include <omp.h>
//#include <pthread.h>

#ifdef WIN32
#include <io.h>
#endif


void init()
{
#ifdef WIN32
	_setmode(_fileno(stdout),O_BINARY);
#endif
	av_register_all();

#ifdef _DEBUG
	av_log_set_level(AV_LOG_INFO);
#else
	av_log_set_level(AV_LOG_FATAL);
#endif

}

char* avpacketMemBuff;


#ifdef WIN32
int round(double val)
{
    return (int)floor(val + 0.5);
}
#endif
int read_av(AVFormatContext *fmt_ctx, 
			AVCodecContext *vc_ctx, 
			AVCodecContext *ac_ctx, 
			AVCodec *vcodec, 
			AVCodec *acodec, 
			int vi, 
			int ai, 
			AVFrame **frame)
{
	AVPacket packet;
	int return_code=0;
	int vframe_finished = 0;
	int aframe_finished = 0;
	int stream_finished = 0;
	int ret;

	while(!vframe_finished && !stream_finished && !aframe_finished)
	{
		av_init_packet(&packet);
		return_code = av_read_frame(fmt_ctx, &packet);//it can be video or audio
		if(return_code<0)
		{
			packet.size = 0;
			packet.data = NULL;

			if(packet.stream_index==vi)
				ret = avcodec_decode_video2(vc_ctx, *frame, &vframe_finished, &packet);
			if(packet.stream_index==ai)
				ret = avcodec_decode_audio4(ac_ctx, *frame, &aframe_finished, &packet);

			if(!vframe_finished)
			{
				stream_finished = 1;
			}

			if(return_code == -EAGAIN)
				stream_finished = 0;

		} 
		else if(packet.stream_index==vi)
		{
			ret = avcodec_decode_video2(vc_ctx, *frame, &vframe_finished, &packet);
		}
		else if(packet.stream_index==ai)
		{
			avcodec_decode_audio4(ac_ctx, *frame, &aframe_finished, &packet);
//			(*frame)->pict_type = vframe->pict_type;
//			(*frame)->pict_type = AV_PICTURE_TYPE_NONE;
//			av_freep(&vframe);
		}

		av_free_packet(&packet);
	}
	return(stream_finished);
}
int read_av_keyFrame(AVFormatContext *fmt_ctx, 
			AVCodecContext *vc_ctx, 
			AVCodecContext *ac_ctx, 
			AVCodec *vcodec, 
			AVCodec *acodec, 
			int vi, 
			int ai, 
			AVFrame **frame)
{
	AVPacket packet;
	int return_code=0;
	int vframe_finished = 0;
	int aframe_finished = 0;
	int stream_finished = 0;
	int ret;

	while(!vframe_finished && !stream_finished && !aframe_finished)
	{
		av_init_packet(&packet);
		return_code = av_read_frame(fmt_ctx, &packet);//it can be video or audio
		if(return_code<0)
		{
			packet.size = 0;
			packet.data = NULL;

			if(packet.stream_index==vi)
				ret = avcodec_decode_video2(vc_ctx, *frame, &vframe_finished, &packet);
			if(packet.stream_index==ai)
				ret = avcodec_decode_audio4(ac_ctx, *frame, &aframe_finished, &packet);

			if(!vframe_finished)
			{
				stream_finished = 1;
			}

			if(return_code == -EAGAIN)
				stream_finished = 0;

		} 
		else if(  (packet.stream_index==vi) && (packet.flags&AV_PKT_FLAG_KEY)!=0)
		{
			ret = avcodec_decode_video2(vc_ctx, *frame, &vframe_finished, &packet);
		}
		else if(packet.stream_index==ai)
		{
			//avcodec_decode_audio4(ac_ctx, *frame, &aframe_finished, &packet);
//			(*frame)->pict_type = vframe->pict_type;
			(*frame)->pict_type = AV_PICTURE_TYPE_NONE;
//			av_freep(&vframe);
		}
		av_free_packet(&packet);
	}
	return(stream_finished);
}
int open_av(char filename[],
                        AVFormatContext *fmt_ctx,
                        AVCodecContext **vc_ctx,
                        AVCodecContext **ac_ctx,
                        AVCodec **vcodec,
                        AVCodec **acodec,
                        int *vi,
                        int *ai,
                        int *width,
                        int *height,
                        double *fps,
                        double *duration
                        )
{
	return open_av_headers(filename,fmt_ctx,vc_ctx,ac_ctx,vcodec,acodec,vi,ai,width,height,fps,duration,"");
}

int open_av_headers(char filename[], 
			AVFormatContext *fmt_ctx, 
			AVCodecContext **vc_ctx, 
			AVCodecContext **ac_ctx, 
			AVCodec **vcodec, 
			AVCodec **acodec, 
			int *vi,			
			int *ai,
			int *width, 
			int *height, 
			double *fps, 
			double *duration,
			char *httpheaders 
			)
{

	int i;
	if(filename=="nothing"){
		if(avformat_open_input(&fmt_ctx, NULL, NULL, NULL) != 0)//args 3,4 are NULL if opening a file
        	{
                	fprintf(stderr, "ERROR: Invalid input file %s\n", filename);
	                return(1);
        	}
	}else{
		if(strcmp(httpheaders,"")==0){
			//if(avformat_open_input(&fmt_ctx, filename, NULL, NULL) != 0)//args 3,4 are NULL if opening a file
        	        //{
                        //	fprintf(stderr, "ERROR: Invalid input file %s\n", filename);
	                //        return(1);
                	//}
                		//av_log(NULL, AV_LOG_WARNING, "begin avformat open input\n");
                        AVDictionary *opts = NULL;
						av_dict_set(&opts,"timeout","500000",0);
						//av_dict_set(&opts,"timeout","200000",0);
						//av_dict_set(&opts,"timeout","500000",0);
                        //av_dict_set(&opts, "timeout", "500000", 0);
                        if(avformat_open_input(&fmt_ctx, filename, NULL, &opts) != 0)//args 3,4 are NULL if opening a file
                        {
                                fprintf(stderr, "ERROR: Invalid input file %s\n", filename);
								av_dict_free(&opts);
                                return(1);
                        }
						//av_log(NULL, AV_LOG_WARNING, "avformat open input over\n");
						av_dict_free(&opts);
		}else{
			char headers[1024] = {0};
			//snprintf(headers, 1024, "Host: %s\r\n",httpheaders);
			snprintf(headers, 1024, "%s\r\n",httpheaders);
			AVDictionary *opts = NULL;
			//av_dict_set(&opts, "timeout", "5000", 0);
			av_dict_set(&opts, "headers", headers, 0);
			
			//fprintf(stderr, "util: %s\nheaders: %s \n",filename,headers);
			if(avformat_open_input(&fmt_ctx, filename, NULL, &opts) != 0)//args 3,4 are NULL if opening a file
                        {
                                fprintf(stderr, "ERROR: Invalid input file %s\n", filename);
                                return(1);
                        }
			av_dict_free( &opts);
		}
	}
	//fmt_ctx->probesize = 10000*1024;
	//fmt_ctx->max_analyze_duration = 5*AV_TIME_BASE;//3M
	fmt_ctx->probesize = 100000*1024;
	fmt_ctx->max_analyze_duration = 10*AV_TIME_BASE;//3M
	if(avformat_find_stream_info(fmt_ctx, NULL) < 0)
	{
		fprintf(stderr, "ERROR: unable to find stream info\n");
		return(1);
	}
	//av_log(NULL, AV_LOG_WARNING, "identifying video stream\n");
/*#ifdef _DEBUG
	av_dump_format(fmt_ctx, 0, filename, 0);
#endif*/
	//av_log(NULL, AV_LOG_WARNING, "video format stream num:", (int)(fmt_ctx->nb_streams));
	for(i=0; i< (int)(fmt_ctx->nb_streams) ; i++)//identifying video stream
	{
		if(fmt_ctx->streams[i]->codec->codec_type==AVMEDIA_TYPE_VIDEO)
		{
			*vi = i;
			*vc_ctx = fmt_ctx->streams[i]->codec;
		}
		else if(fmt_ctx->streams[i]->codec->codec_type==AVMEDIA_TYPE_AUDIO)
		{
			*ai = i;
			*ac_ctx = fmt_ctx->streams[i]->codec;
		}
	}
	//av_log(NULL, AV_LOG_WARNING, "begin find decoder\n");
	if(*vc_ctx!=NULL)
		*vcodec = avcodec_find_decoder((*vc_ctx)->codec_id);
#ifdef _DEBUG
	else
	{
		fprintf(stderr, "WARNING: Video stream not found\n");
//		return(1);
	}
#endif

	if(*ac_ctx!=NULL)
		*acodec = avcodec_find_decoder((*ac_ctx)->codec_id);
#ifdef _DEBUG
	else
	{
		fprintf(stderr, "WARNING: Audio stream not found\n");
//		return(1);
	}
#endif
	//av_log(NULL, AV_LOG_WARNING, "find decoder over\n");
	*duration = (double)fmt_ctx->duration/AV_TIME_BASE;
	if(*duration < 0)
		*duration = 0;
	if(*vi>=0)
	{
		//*fps = (double)((*vc_ctx)->time_base.den)/((*vc_ctx)->time_base.num)/(*vc_ctx)->ticks_per_frame;
		// return average fps
		*fps = (double)(fmt_ctx->streams[*vi]->avg_frame_rate.num)/fmt_ctx->streams[*vi]->avg_frame_rate.den;
		*height = (*vc_ctx)->height;
		*width = (*vc_ctx)->width;
	}

#ifdef _DEBUG
	if(*vi>=0)
	{
		fprintf(stderr,"\nVideo codec name: %s\n", (*vcodec)->name);
		fprintf(stderr,"Video framerate (header): %3.2f fps (%d, %d, %d)\n",*fps, (*vc_ctx)->time_base.den, (*vc_ctx)->time_base.num, (*vc_ctx)->ticks_per_frame);
		fprintf(stderr,"Video average framerate : %3.4f (%d / %d)\n", (double)(fmt_ctx->streams[*vi]->avg_frame_rate.num)/fmt_ctx->streams[*vi]->avg_frame_rate.den, fmt_ctx->streams[*vi]->avg_frame_rate.num, fmt_ctx->streams[*vi]->avg_frame_rate.den);
		fprintf(stderr,"Video resolution: %dx%d\n",*width,*height );
		fprintf(stderr,"Duration = %3.2f\n",*duration);
	}
	if(*ai>=0)
	{
		fprintf(stderr,"\nAudio codec name: %s\n", (*acodec)->name);
		fprintf(stderr,"Audio sample rate: %dHz\n", (*ac_ctx)->sample_rate);
		fprintf(stderr,"Audio channels: %d\n", (*ac_ctx)->channels);
		fprintf(stderr,"Audio channel layout: %d\n", (*ac_ctx)->channel_layout);
	}

#endif
	if(*vc_ctx!=NULL)
	{
		//(*vc_ctx)->thread_count = 1;
		(*vc_ctx)->thread_count = getNumber_of_thread();
//printf("Util Number_thread:%d\n",getNumber_of_thread());
		//av_log(NULL, AV_LOG_WARNING, "begin avcodec open video\n");
		if(avcodec_open2(*vc_ctx, *vcodec, NULL)<0)
		{
			fprintf(stderr, "ERROR: unsupported video codec\n");
			return(19);
		}
		//av_log(NULL, AV_LOG_WARNING, "avcodec open video over\n");
	}
	if(*ac_ctx!=NULL){
		if(avcodec_open2(*ac_ctx, *acodec, NULL)<0)
		{
			fprintf(stderr, "WARNING: unsupported audio codec\n");
		}
	}
	return(0);
}

int set_avfilter(AVFilterGraph *filter_graph, AVFormatContext *fmt_ctx, AVStream *vstream, AVStream *astream, int sampleRate)
{
	static char strbuf[512];
	static AVFilterContext *abuffer_ctx = NULL;
	static AVFilterContext *volume_ctx = NULL;
	static AVFilterContext *aformat_ctx = NULL;
	static AVFilterContext *abuffersink_ctx = NULL;

	AVFilter *abuffer;
	//AVFilter *volume;
	AVFilter *aformat;
	AVFilter *abuffersink;
    AVCodecContext *avctx = astream->codec;
    AVRational time_base = astream->time_base;
	int err;
	uint64_t channel_layout;
	uint64_t out_chlayout;
	avfilter_register_all();
	abuffer = avfilter_get_by_name("abuffer");
	//volume = avfilter_get_by_name("volumedetect");
	aformat = avfilter_get_by_name("aformat");
	abuffersink = avfilter_get_by_name("abuffersink");
	
	if(avctx->channel_layout==0)
		channel_layout = av_get_default_channel_layout(avctx->channels);
	else
		channel_layout = avctx->channel_layout;
// create abuffer filter
#ifdef WIN32
	sprintf_s(strbuf, sizeof(strbuf), "time_base=%d/%d:sample_rate=%d:sample_fmt=%s:channel_layout=0x%"PRIx64, time_base.num, time_base.den, avctx->sample_rate, av_get_sample_fmt_name(avctx->sample_fmt), channel_layout);
#else
	snprintf(strbuf, sizeof(strbuf), "time_base=%d/%d:sample_rate=%d:sample_fmt=%s:channel_layout=0x%"PRIx64, time_base.num, time_base.den, avctx->sample_rate, av_get_sample_fmt_name(avctx->sample_fmt), channel_layout);
#endif	
    err = avfilter_graph_create_filter(&abuffer_ctx, abuffer, "in", strbuf, NULL, filter_graph);
	if (err < 0) {
        av_log(NULL, AV_LOG_ERROR, "error initializing abuffer filter\n");
        return err;
    }

// create volume filter
//	sprintf_s(strbuf, sizeof(strbuf),"volume=-20dB");
    //err = avfilter_graph_create_filter(&volume_ctx, volume, "volumedetect", NULL/*strbuf*/, NULL, filter_graph);
    //if (err < 0) {
    //    av_log(NULL, AV_LOG_ERROR, "error initializing volume filter\n");
    //    return err;
    //}
// create aformat filter
	out_chlayout = AV_CH_LAYOUT_STEREO;
#ifdef WIN32
    //sprintf_s(strbuf, sizeof(strbuf), "sample_fmts=%s:sample_rates=%d:channel_layouts=0x%"PRIx64, av_get_sample_fmt_name(AV_SAMPLE_FMT_FLT), 8000, out_chlayout);
	sprintf_s(strbuf, sizeof(strbuf), "sample_fmts=%s:sample_rates=%d:channel_layouts=0x%"PRIx64, av_get_sample_fmt_name(AV_SAMPLE_FMT_FLT), sampleRate, channel_layout);
#else
    snprintf(strbuf, sizeof(strbuf), "sample_fmts=%s:sample_rates=%d:channel_layouts=0x%"PRIx64, av_get_sample_fmt_name(AV_SAMPLE_FMT_FLT), sampleRate, channel_layout);
#endif
//    fprintf(stderr, "aformat: %s\n", strbuf);
    err = avfilter_graph_create_filter(&aformat_ctx, aformat,
            NULL, strbuf, NULL, filter_graph);
    if (err < 0) {
        av_log(NULL, AV_LOG_ERROR, "unable to create aformat filter\n");
        return err;
    }

 // create abuffersink filter
    err = avfilter_graph_create_filter(&abuffersink_ctx, abuffersink, "out", NULL, NULL, filter_graph);
    if (err < 0) {
        av_log(NULL, AV_LOG_ERROR, "unable to create aformat filter\n");
        return err;
    }
    if(err >= 0) err = avfilter_link(abuffer_ctx, 0, aformat_ctx, 0);
    //if (err >= 0) err = avfilter_link(volume_ctx, 0, aformat_ctx, 0);
    if (err >= 0) err = avfilter_link(aformat_ctx, 0, abuffersink_ctx, 0);
    if (err < 0) {
        av_log(NULL, AV_LOG_ERROR, "error connecting filters\n");
        return err;
    }

    err = avfilter_graph_config(filter_graph, NULL);
    if (err < 0) {
        av_log(NULL, AV_LOG_ERROR, "error configuring the filter graph\n");
        return err;
    }
    return 0;
}

unsigned char **allocate_2d_uchar(int y, int x)
{
	int j;
	unsigned char **mem_blk;

	mem_blk = malloc(sizeof(unsigned char*)*y);
	if(mem_blk==NULL){
		return(NULL);
	}
	for(j=0; j<y; j++)
	{
		mem_blk[j] = malloc(sizeof(unsigned char)*x);
		if(mem_blk[j]==NULL){
			int i;
			for (i=0; i< j; i++){
				free(mem_blk[i]);
			}
			free(mem_blk);
			return(NULL);
			break;
		}
	}
	return(mem_blk);
}

double **allocate_2d_double(int y, int x)
{
	int j;
	double **mem_blk;

	mem_blk = malloc(sizeof(double*)*y);
	if(mem_blk==NULL)
		return(NULL);
	for(j=0; j<y; j++)
	{
		mem_blk[j] = malloc(sizeof(double)*x);
		if(mem_blk[j]==NULL){
			int i;
			for (i=0; i< j; i++){
				free(mem_blk[i]);
			}
			free(mem_blk);
			return(NULL);
		}
	}
	return(mem_blk);
}
float **allocate_2d_float(int y, int x)
{
	int j;
	float **mem_blk;

	mem_blk = malloc(sizeof(float*)*y);
	if(mem_blk==NULL)
		return(NULL);
	for(j=0; j<y; j++)
	{
		mem_blk[j] = malloc(sizeof(float)*x);
		if(mem_blk[j]==NULL){
			int i;
			for (i=0; i< j; i++){
				free(mem_blk[i]);
			}
			free(mem_blk);
			return(NULL);
		}
	}
	return(mem_blk);
}
int **allocate_2d_int(int y, int x)
{
	int j;
	int **mem_blk;

	mem_blk = malloc(sizeof(int*)*y);
	if(mem_blk==NULL)
		return(NULL);
	for(j=0; j<y; j++)
	{
		mem_blk[j] = malloc(sizeof(int)*x);
		if(mem_blk[j]==NULL){
			int i;
			for (i=0; i< j; i++){
				free(mem_blk[i]);
			}
			free(mem_blk);
			return(NULL);
		}
	}
	return(mem_blk);
}

unsigned char*** allocate_3d_uchar(int t, int x, int y)
{
	unsigned char*** output3d;
	int i, j;
	output3d = malloc(t*sizeof(unsigned char**));
	if (output3d==NULL){
		return(NULL);
	}
	//if(output3d==NULL) { fprintf(stderr,"Out of memory"); return(NULL);}
	for(j = 0; j<t; j++)
	{
		output3d[j] = malloc(x*sizeof(unsigned char*));
		if (output3d[j]==NULL){
			int s,t;
			for (s=0; s< j; s++){
				for (t=0; t<x; t++){
					free(output3d[s][t]);
				}
				free(output3d[s]);
			}
			free(output3d);
			return(NULL);
		}
		//if(output3d[j]==NULL) { fprintf(stderr,"Out of memory"); return(NULL);}
		for(i = 0; i<x; i++)
		{
			output3d[j][i] = malloc(y*sizeof(unsigned char));
			if (output3d[j][i]==NULL){
				int s;
				int t;
				for (s=0; s< j; s++){
					for (t=0; t<x; t++){
						free(output3d[s][t]);
					}
					free(output3d[s]);
				}
				for (t=0; t< i; t++){
					free(output3d[s][t]);
				}
				free(output3d[j]);
				free(output3d);
				return(NULL);
			}
			//if(output3d[j][i]==NULL) { fprintf(stderr,"Out of memory"); return(NULL);}
		}
	}
	return(output3d);
}
double sqr(double x)
{
	return(x*x);
}

double in_the_ellipse(int y, int x, int cy, int cx, double sy, double sx, double theta)
{
	const double pi=3.14159;
	double xx, yy, cxx, cyy;
	
	xx = x * cos(theta/180*pi) - y * sin(theta/180*pi);//rotation
	yy = x * sin(theta/180*pi) + y * cos(theta/180*pi);
	cxx = cx * cos(theta/180*pi) - cy * sin(theta/180*pi);
	cyy = cx * sin(theta/180*pi) + cy * cos(theta/180*pi);
	
	if( sqr((cxx-xx)/sx) + sqr((cyy-yy)/sy) < 1 )
		return( sqr(sqr((cxx-xx)/sx) + sqr((cyy-yy)/sy)) );
	else return(1.0);
}


int yuv_image_adjust(AVFrame *frame, 
					 int auto_gamma, 
					 double brightness, //0.0 = no change
					 double contrast,	//1.0 = no change
					 double satuation)	//1.0 = no change
{

	int i, j;
	int u, v;
	int mapY[256];
	int histY[256];
	int mapU[256][256];
	int mapV[256][256];
	double sat_factor;
	double gamma;
	int acc;
	int median, mean;

	//Y
	for(i=0; i<256; i++)
		histY[i] = 0;

	for(j=0; j<frame->height; j++)
		for(i=0; i<frame->width; i++)
			histY[*(frame->data[0]+(frame->linesize[0]*j)+i)]++;

	median=0;
	acc=0;
	while(acc < frame->height*frame->width/2)// median of Y component
		acc += histY[median++];

	gamma = 1.0;
	if(auto_gamma)
	{
		if(median < 100)
		{
			gamma = median / 100;
			if(gamma<0.75)
				gamma = 0.75;
		}
		else if(median > 180)
		{
			gamma = (median-80)/100;
			if(gamma>1.25)
				gamma = 1.25;
		}
		for(i=0; i<256; i++)
			mapY[i] = (int)(255 * pow(i/256.0, gamma));
	} 
	else
	{
		for(i=0; i<256; i++)
			mapY[i] = i;
	}


	acc=0;
	for(i=0; i<256; i++)
		acc += histY[i]*mapY[i];
	mean = acc / frame->height / frame->width;

	for(i=0; i<256; i++)
	{
		mapY[i] = (int)(((mapY[i]-mean) * contrast ) + brightness + mean);
		if(mapY[i]>255)
			mapY[i]=255;
		else if(mapY[i]<0)
			mapY[i]=0;
	}

	for(j=0; j<frame->height; j++)//contrast, brightness
		for(i=0; i<frame->width; i++)
			*(frame->data[0]+(frame->linesize[0]*j)+i) = (uint8_t)mapY[*(frame->data[0]+(frame->linesize[0]*j)+i)];

	//UV
	for(u=0; u<256; u++)//U
		for(v=0; v<256; v++)//V
		{

			if(u>125 && u<131 && v>125 && v<131)//almost white
			{
				mapU[u][v] = u;
				mapV[u][v] = v;
			}
			else if( in_the_ellipse(v, u, 146, 108, 11.0+8.0, 32.0+8.0, 45.0) < 1.0)//extended region for more natural tuning
			{
//				mapU[u][v] = (int)((u-128)*sqrt(sqrt(satuation))+128);
//				mapV[u][v] = (int)((v-128)*sqrt(sqrt(satuation))+128);
//				mapU[u][v] = 0;
//				mapV[u][v] = 0;
				sat_factor = in_the_ellipse(v, u, 146, 108, 11.0+8.0, 32.0+8.0, 45.0);
				mapU[u][v] = (int)((u-128)*( (satuation-1)*sat_factor+1 )+128);
				mapV[u][v] = (int)((v-128)*( (satuation-1)*sat_factor+1 )+128);
			} 
			else
			{
				mapU[u][v] = (int)((u-128)*satuation+128);
				mapV[u][v] = (int)((v-128)*satuation+128);

				if(mapU[u][v]>255)
					mapU[u][v]=255;
				else if(mapU[u][v]<0)
					mapU[u][v]=0;
				if(mapV[u][v]>255)
					mapV[u][v]=255;
				else if(mapV[u][v]<0)
					mapV[u][v]=0;
			}
		}
	//UV
	for(j=0; j<frame->height/2; j++)//U' = ((U - 128) x C x S) + 128
		for(i=0; i<frame->width/2; i++)
		{
			u = *(frame->data[1]+(frame->linesize[1]*j)+i);
			v = *(frame->data[2]+(frame->linesize[2]*j)+i);

			*(frame->data[1]+(frame->linesize[1]*j)+i) = (uint8_t)mapU[u][v];
			*(frame->data[2]+(frame->linesize[2]*j)+i) = (uint8_t)mapV[u][v];
		}

	return(0);
}

int to_image(AVFrame* frame, char filename[])
{
	FILE *f;
	AVCodec *codec;
	AVCodecContext *codec_ctx=NULL;
	AVPacket packet;
	int got_packet;	

	codec = avcodec_find_encoder(AV_CODEC_ID_MJPEG);
	if(!codec)
	{
		fprintf(stderr, "No JPEG encoder in this build\n");
		exit(1);
	}

	codec_ctx = avcodec_alloc_context3(codec);
	codec_ctx->qmax=5;
	codec_ctx->qmin=5;
	codec_ctx->width = frame->width;
	codec_ctx->height = frame->height;
	codec_ctx->pix_fmt = AV_PIX_FMT_YUVJ420P;
	codec_ctx->time_base.den = 1;
	codec_ctx->time_base.num = 1;
	avcodec_open2(codec_ctx, codec, NULL);
#ifdef WIN32
	fopen_s(&f, filename, "wb");
#else
	f = fopen(filename, "wb");
#endif

	av_init_packet(&packet);
	packet.data = NULL;
	packet.size = 0;

	avcodec_encode_video2(codec_ctx, &packet, frame, &got_packet);
	fwrite(packet.data, 1, packet.size, f);
	av_free_packet(&packet);
	avcodec_close(codec_ctx);
	fflush(f);
	fclose(f);
	return(0);
}

int to_image_png(AVFrame* frame, char filename[])
{
	FILE *f;
	AVCodec *codec;
	AVCodecContext *codec_ctx=NULL;
	AVPacket packet;
	int got_packet;	

	codec = avcodec_find_encoder(AV_CODEC_ID_PNG);
	if(!codec)
	{
		fprintf(stderr, "No PNG encoder in this build of libavcodec\n");
		exit(1);
	}

	codec_ctx = avcodec_alloc_context3(codec);
	codec_ctx->qmax=5;
	codec_ctx->qmin=5;
	codec_ctx->width = frame->width;
	codec_ctx->height = frame->height;
	codec_ctx->pix_fmt = frame->format;
	codec_ctx->time_base.den = 1;
	codec_ctx->time_base.num = 1;
	avcodec_open2(codec_ctx, codec, NULL);
#ifdef WIN32
	fopen_s(&f, filename, "wb");
#else
	f = fopen(filename, "wb");
#endif

	av_init_packet(&packet);
	packet.data = NULL;
	packet.size = 0;

	avcodec_encode_video2(codec_ctx, &packet, frame, &got_packet);
	fwrite(packet.data, 1, packet.size, f);
	av_free_packet(&packet);
	
	avcodec_close(codec_ctx);
	
	fflush(f);
	fclose(f);
	return(0);
}

int to_image_2(AVFrame* frame, 
			   char filename[], 
			   int width, 
			   int height, 
			   int crop, 
			   int encoder, 
			   int auto_gamma,
			   double brightness, 
			   double contrast, 
			   double satuation)//crop, pad, resize, adjust and encode image
{
	FILE *f;
	AVCodec *codec;
	AVCodecContext *codec_ctx=NULL;
	AVPacket packet;
	int got_packet;
	int i;

	int cp_top, cp_bottom, cp_left, cp_right;
	AVFrame *cp_frame;
	uint8_t *cp_buffer;
	int cp_buffer_size;

	AVFrame *zoomed_frame;
	uint8_t *zoomed_buffer;
	int zoomed_buffer_size;
	int pix_fmt;

	struct SwsContext *sws_ctx;
	struct SwsFilter *sws_filter;

	if(encoder!=AV_CODEC_ID_MJPEG && encoder!=AV_CODEC_ID_PNG)
	{
		fprintf(stderr, "Unsupported image format codec_id=%d",encoder);
		return(1);
	}

	cp_top = 0;
	cp_bottom = 0;
	cp_left = 0;
	cp_right = 0;
	cp_frame = av_frame_alloc();
	if(crop)
	{
		if( frame->height*width > frame->width*height )//top bottom cropping
		{
			cp_top = ( frame->height - (frame->width * height / width) ) / 2;
			cp_bottom = frame->height - cp_top;
			cp_top += cp_top%2;
			cp_bottom += cp_bottom%2;
			cp_left = 0;
			cp_right = frame->width;
		} 
		else
		{
			cp_left = ( frame->width - (frame->height * width / height) ) / 2;
			cp_right = frame->width - cp_left;
			cp_left += cp_left%2;
			cp_right += cp_right%2;
			cp_top = 0;
			cp_bottom = frame->height;
		}
		cp_buffer_size = avpicture_get_size(AV_PIX_FMT_YUV420P, cp_right-cp_left, cp_bottom-cp_top );
		cp_buffer = (uint8_t *)av_malloc(cp_buffer_size*sizeof(uint8_t));
		avpicture_fill((AVPicture *)(cp_frame), cp_buffer, AV_PIX_FMT_YUV420P, cp_right-cp_left, cp_bottom-cp_top );
		cp_frame->width = cp_right-cp_left;
		cp_frame->height = cp_bottom-cp_top;

		for(i=0; i<cp_frame->height; i++)
		{
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*i, frame->data[0]+frame->linesize[0]*(i+cp_top)+cp_left, cp_frame->width);
		}
		for(i=0; i<cp_frame->height/2; i++)
		{
			memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, frame->data[1]+frame->linesize[1]*(i+cp_top/2)+cp_left/2, cp_frame->width/2);
			memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, frame->data[2]+frame->linesize[1]*(i+cp_top/2)+cp_left/2, cp_frame->width/2);
		}
	} else //pad
	{
		if( frame->height*width > frame->width*height )//left right padding
		{
			cp_left = ((frame->height * width / height) - frame->width)/2;
			cp_right = (frame->height * width / height);
			cp_left += cp_left%2;
			cp_right += cp_right%2;
			cp_top = 0;
			cp_bottom = frame->height;
		} 
		else
		{
			cp_top = ((frame->width * height / width) - frame->height)/2;
			cp_bottom = (frame->width * height / width);
			cp_top += cp_top%2;
			cp_bottom += cp_bottom%2;
			cp_left = 0;
			cp_right = frame->width;
		}
		cp_buffer_size = avpicture_get_size(AV_PIX_FMT_YUV420P, cp_right, cp_bottom );
		cp_buffer = (uint8_t *)av_malloc(cp_buffer_size*sizeof(uint8_t));
		avpicture_fill((AVPicture *)(cp_frame), cp_buffer, AV_PIX_FMT_YUV420P, cp_right, cp_bottom );
		cp_frame->width = cp_right;
		cp_frame->height = cp_bottom;

		for(i=0; i<cp_frame->width; i++)
			*(cp_frame->data[0]+i)=0;
		for(i=0; i<cp_frame->height; i++)
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*i, cp_frame->data[0], cp_frame->width);
		for(i=0; i<cp_frame->width/2; i++)
		{
			*(cp_frame->data[1]+i)=128;
			*(cp_frame->data[2]+i)=128;
		}
		for(i=0; i<cp_frame->height/2; i++)
		{
			memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, cp_frame->data[1], cp_frame->width/2);
			memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, cp_frame->data[1], cp_frame->width/2);
		}
		for(i=0; i<frame->height; i++)
		{
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*(i+cp_top)+cp_left, frame->data[0]+frame->linesize[0]*i, frame->width);
		}
		for(i=0; i<frame->height/2; i++)
		{
			memcpy(cp_frame->data[1]+cp_frame->linesize[1]*(i+cp_top/2)+cp_left/2, frame->data[1]+frame->linesize[1]*i, frame->width/2);
			memcpy(cp_frame->data[2]+cp_frame->linesize[2]*(i+cp_top/2)+cp_left/2, frame->data[2]+frame->linesize[2]*i, frame->width/2);
		}
	}//done cropping 
	
	yuv_image_adjust(cp_frame, auto_gamma, brightness, contrast, satuation);

	//start zooming
	zoomed_frame = av_frame_alloc();
	pix_fmt = (encoder==AV_CODEC_ID_MJPEG)?AV_PIX_FMT_YUVJ420P:AV_PIX_FMT_RGB24;
	zoomed_buffer_size = avpicture_get_size( pix_fmt, width, height );
	zoomed_buffer = (uint8_t *)av_malloc(zoomed_buffer_size*sizeof(uint8_t));
	avpicture_fill((AVPicture *)(zoomed_frame), zoomed_buffer, pix_fmt, width, height );
	zoomed_frame->width = width;
	zoomed_frame->height = height;
	zoomed_frame->format = pix_fmt;
	sws_filter = sws_getDefaultFilter(0,0,1,1,0,0,1);
	sws_filter->lumH = sws_getConstVec(1,3);
	sws_filter->lumH->coeff[0] = -1;
	sws_filter->lumH->coeff[1] = 3;
	sws_filter->lumH->coeff[2] = -1;

	sws_filter->lumV = sws_getConstVec(1,3); 
	sws_filter->lumV->coeff[0] = -1;
	sws_filter->lumV->coeff[1] = 3;
	sws_filter->lumV->coeff[2] = -1;


	sws_filter->chrH = sws_getConstVec(1,3); 
	sws_filter->chrV = sws_getConstVec(1,3); 
	

	sws_ctx = sws_getContext(cp_frame->width, cp_frame->height, AV_PIX_FMT_YUV420P, width, height,pix_fmt, SWS_SPLINE, NULL, NULL, NULL);

	if(sws_scale(sws_ctx, (const uint8_t **)(cp_frame->data), cp_frame->linesize, 0, cp_frame->height, zoomed_frame->data, zoomed_frame->linesize)!=zoomed_frame->height)
	{
		fprintf(stderr, "ERROR in SWS\n");
		return(1);
	}
	//done zooming
	codec = avcodec_find_encoder(encoder);
	if(!codec)
	{
		fprintf(stderr, "No %s encoder in this build\n", (encoder==AV_CODEC_ID_MJPEG)?"JPEG":"PNG");
		exit(1);
	}

	codec_ctx = avcodec_alloc_context3(codec);
	codec_ctx->qmax=2;
	codec_ctx->qmin=2;

	codec_ctx->width = zoomed_frame->width;
	codec_ctx->height = zoomed_frame->height;
	codec_ctx->pix_fmt = pix_fmt;
	codec_ctx->time_base.den = 1;
	codec_ctx->time_base.num = 1;
	avcodec_open2(codec_ctx, codec,NULL);

#ifdef WIN32
	fopen_s(&f, filename, "wb");
#else
	f = fopen(filename, "wb");
#endif
	av_init_packet(&packet);
	packet.data = NULL;
	packet.size = 0;

	avcodec_encode_video2(codec_ctx, &packet, zoomed_frame, &got_packet);
	fwrite(packet.data, 1, packet.size, f);
	av_free_packet(&packet);
	avcodec_close(codec_ctx);
	fflush(f);
	fclose(f);
	return(0);
}

int to_image_3(AVFrame* frame, 
			   char filename[], 
			   int width, 
			   int height, 
			   int crop, 
			   int encoder,
			   int auto_gamma,
			   double brightness, 
			   double contrast, 
			   double satuation,
			   AVPacket *img_pkt,
			   int crf)//crop, pad, resize, adjust and encode image
{
//	av_log_set_level(AV_LOG_QUIET);

	AVFormatContext *format_ctx=NULL;
//	AVOutputFormat *format=NULL;
	AVCodecContext *codec_ctx=NULL;
	AVCodec *codec=NULL;
	AVPacket packet;
	AVDictionary *opts=NULL;
	AVStream *out_stream;
	char optstr[1024];

	int got_packet;
	int i;

	int cp_top, cp_bottom, cp_left, cp_right;
	AVFrame *cp_frame;
	uint8_t *cp_buffer;
	int cp_buffer_size;

	AVFrame *zoomed_frame;
	uint8_t *zoomed_buffer;
	int zoomed_buffer_size;
	int pix_fmt;

	struct SwsContext *sws_ctx;
	int write_flag = 1;

	if(crf<1)
		crf=28;

	if(encoder!=AV_CODEC_ID_MJPEG && encoder!=AV_CODEC_ID_PNG && encoder!=AV_CODEC_ID_HEVC && encoder!=AV_CODEC_ID_H264 && encoder!=AV_CODEC_ID_VP9)
	{
		fprintf(stderr, "Unsupported image format codec_id=%d\n",encoder);
		return(1);
	}

	format_ctx = avformat_alloc_context();
	format_ctx->oformat = av_guess_format(NULL, filename, NULL);
	if(format_ctx->oformat == 0)
	{
		//fprintf(stderr, "Unsupported output format %s, not writing\n", filename);
		write_flag = 0;
	}
	cp_top = 0;
	cp_bottom = 0;
	cp_left = 0;
	cp_right = 0;
	cp_frame = av_frame_alloc();

	if(crop)
	{
		if( frame->height*width > frame->width*height )//top bottom cropping
		{
			cp_top = ( frame->height - (frame->width * height / width) ) / 2;
			cp_bottom = frame->height - cp_top;
			cp_top += cp_top%2;
			cp_bottom += cp_bottom%2;
			cp_left = 0;
			cp_right = frame->width;
		} 
		else
		{
			cp_left = ( frame->width - (frame->height * width / height) ) / 2;
			cp_right = frame->width - cp_left;
			cp_left += cp_left%2;
			cp_right += cp_right%2;
			cp_top = 0;
			cp_bottom = frame->height;
		}
		cp_frame->width = cp_right-cp_left;
		cp_frame->height = cp_bottom-cp_top;
		cp_frame->format = frame->format;
		cp_buffer_size = avpicture_get_size(frame->format, cp_right-cp_left, cp_bottom-cp_top );
		cp_buffer = (uint8_t *)av_malloc(cp_buffer_size*sizeof(uint8_t));

		avpicture_fill((AVPicture *)(cp_frame), cp_buffer, frame->format, cp_right-cp_left, cp_bottom-cp_top );

		for(i=0; i<cp_frame->height; i++)
		{
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*i, frame->data[0]+frame->linesize[0]*(i+cp_top)+cp_left, cp_frame->width);
		}
		if(frame->format==AV_PIX_FMT_YUV420P || frame->format==AV_PIX_FMT_YUVJ420P)
		{
			for(i=0; i<cp_frame->height/2; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, frame->data[1]+frame->linesize[1]*(i+cp_top/2)+cp_left/2, cp_frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, frame->data[2]+frame->linesize[1]*(i+cp_top/2)+cp_left/2, cp_frame->width/2);
			}
		}
		else if(frame->format==AV_PIX_FMT_YUV422P || frame->format==AV_PIX_FMT_YUVJ422P)
		{
			for(i=0; i<cp_frame->height; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, frame->data[1]+frame->linesize[1]*(i+cp_top)+cp_left/2, cp_frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, frame->data[2]+frame->linesize[2]*(i+cp_top)+cp_left/2, cp_frame->width/2);
			}
		}
	} else //pad
	{
		if( frame->height*width > frame->width*height )//left right padding
		{
			cp_left = ((frame->height * width / height) - frame->width)/2;
			cp_right = (frame->height * width / height);
			cp_left += cp_left%2;
			cp_right += cp_right%2;
			cp_top = 0;
			cp_bottom = frame->height;
		} 
		else
		{
			cp_top = ((frame->width * height / width) - frame->height)/2;
			cp_bottom = (frame->width * height / width);
			cp_top += cp_top%2;
			cp_bottom += cp_bottom%2;
			cp_left = 0;
			cp_right = frame->width;
		}
		cp_buffer_size = avpicture_get_size(frame->format, cp_right, cp_bottom );
		cp_buffer = (uint8_t *)av_malloc(cp_buffer_size*sizeof(uint8_t));
		avpicture_fill((AVPicture *)(cp_frame), cp_buffer, frame->format, cp_right, cp_bottom );
		cp_frame->width = cp_right;
		cp_frame->height = cp_bottom;

		for(i=0; i<cp_frame->width; i++)//first row
			*(cp_frame->data[0]+i)=0;
		for(i=0; i<cp_frame->height; i++)
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*i, cp_frame->data[0], cp_frame->width);

		if(frame->format==AV_PIX_FMT_YUV420P || frame->format==AV_PIX_FMT_YUVJ420P)
		{
			for(i=0; i<cp_frame->width/2; i++)//first row
			{
				*(cp_frame->data[1]+i)=128;
				*(cp_frame->data[2]+i)=128;
			}
			for(i=0; i<cp_frame->height/2; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, cp_frame->data[1], cp_frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, cp_frame->data[1], cp_frame->width/2);
			}
		}
		else if(frame->format==AV_PIX_FMT_YUV422P || frame->format==AV_PIX_FMT_YUVJ422P)
		{
			for(i=0; i<cp_frame->width/2; i++)//first row
			{
				*(cp_frame->data[1]+i)=128;
				*(cp_frame->data[2]+i)=128;
			}
			for(i=0; i<cp_frame->height; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*i, cp_frame->data[1], cp_frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*i, cp_frame->data[1], cp_frame->width/2);
			}
		}

		for(i=0; i<frame->height; i++)
		{
			memcpy(cp_frame->data[0]+cp_frame->linesize[0]*(i+cp_top)+cp_left, frame->data[0]+frame->linesize[0]*i, frame->width);
		}

		if(frame->format==AV_PIX_FMT_YUV420P || frame->format==AV_PIX_FMT_YUVJ420P)
		{
			for(i=0; i<frame->height/2; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*(i+cp_top/2)+cp_left/2, frame->data[1]+frame->linesize[1]*i, frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*(i+cp_top/2)+cp_left/2, frame->data[2]+frame->linesize[2]*i, frame->width/2);
			}
		}
		else if(frame->format==AV_PIX_FMT_YUV422P || frame->format==AV_PIX_FMT_YUVJ422P)
		{
			for(i=0; i<frame->height; i++)
			{
				memcpy(cp_frame->data[1]+cp_frame->linesize[1]*(i+cp_top)+cp_left/2, frame->data[1]+frame->linesize[1]*i, frame->width/2);
				memcpy(cp_frame->data[2]+cp_frame->linesize[2]*(i+cp_top)+cp_left/2, frame->data[2]+frame->linesize[2]*i, frame->width/2);
			}
		}
			
	}//done cropping 
	
//	yuv_image_adjust(cp_frame, auto_gamma, brightness, contrast, satuation);

	//start zooming
	zoomed_frame = av_frame_alloc();
	if(encoder==AV_CODEC_ID_MJPEG)
		pix_fmt = AV_PIX_FMT_YUVJ420P;
	else if(encoder==AV_CODEC_ID_HEVC || encoder==AV_CODEC_ID_H264 || encoder==AV_CODEC_ID_VP9)
		pix_fmt = AV_PIX_FMT_YUV420P;
	else
		pix_fmt = AV_PIX_FMT_RGB24;

	zoomed_frame->width = width;
	zoomed_frame->height = height;
	zoomed_frame->format = pix_fmt;
	zoomed_buffer_size = avpicture_get_size( pix_fmt, width, height );
	zoomed_buffer = (uint8_t *)av_malloc(zoomed_buffer_size*sizeof(uint8_t));
	avpicture_fill((AVPicture *)(zoomed_frame), zoomed_buffer, pix_fmt, width, height );

	sws_ctx = sws_getContext(cp_frame->width, cp_frame->height, frame->format, width, height,pix_fmt, SWS_SPLINE, NULL, NULL, NULL);

	if(sws_scale(sws_ctx, (const uint8_t **)(cp_frame->data), cp_frame->linesize, 0, cp_frame->height, zoomed_frame->data, zoomed_frame->linesize)!=zoomed_frame->height)
	{
		fprintf(stderr, "ERROR in SWS\n");
		return(1);
	}


	//done zooming
	codec = avcodec_find_encoder(encoder);
	if(!codec)
	{
		if(encoder==AV_CODEC_ID_MJPEG)
			fprintf(stderr, "No JPEG encoder in this build\n");
		if(encoder==AV_CODEC_ID_PNG)
			fprintf(stderr, "No PNG encoder in this build\n");
		if(encoder==AV_CODEC_ID_HEVC)
			fprintf(stderr, "No HEVC encoder in this build\n");
		if(encoder==AV_CODEC_ID_H264)
			fprintf(stderr, "No H264 encoder in this build\n");
		if(encoder==AV_CODEC_ID_VP9)
			fprintf(stderr, "No VP9 encoder in this build\n");
		exit(1);
	}

	codec_ctx = avcodec_alloc_context3(codec);
	if(write_flag)
	{
		out_stream = avformat_new_stream(format_ctx, codec);
		out_stream->codec = codec_ctx;
	}

	codec_ctx->width = zoomed_frame->width;
	codec_ctx->height = zoomed_frame->height;
	codec_ctx->pix_fmt = pix_fmt;
	codec_ctx->time_base.den = 25;
	codec_ctx->time_base.num = 1;

	if(encoder==AV_CODEC_ID_MJPEG)
	{
		codec_ctx->qmax=2;
		codec_ctx->qmin=2;
	}
	else if(encoder==AV_CODEC_ID_HEVC)
	{
	#if (!defined(WIN32))
		snprintf(optstr, 1024, "preset=slow:crf=%d:rect=1:max-merge=3:psy-rdoq=1.00:rdoq-level=2:rd=4:log-level=1",crf);
	#endif	
		av_dict_set(&opts, "x265-params", optstr, 0);	
	}
	else if(encoder==AV_CODEC_ID_H264)
	{
	#if (!defined(WIN32))
		snprintf(optstr, 1024, "crf=%d",crf);
	#endif
		av_dict_set(&opts, "x264opts", optstr, 0);
	}

	avcodec_open2(codec_ctx, codec,&opts);

	if(write_flag)
	{
		avio_open(&format_ctx->pb, filename, AVIO_FLAG_WRITE);
		int retc = avformat_write_header(format_ctx, NULL);
	}
	av_init_packet(&packet);
	packet.data = NULL;
	packet.size = 0;

	avcodec_encode_video2(codec_ctx, &packet, zoomed_frame, &got_packet);
	while(!got_packet)
		avcodec_encode_video2(codec_ctx, &packet, NULL, &got_packet);

	img_pkt->size = packet.size;
	img_pkt->data = av_malloc(packet.size);
	memcpy(img_pkt->data, packet.data, packet.size);

	if(write_flag)
		av_interleaved_write_frame(format_ctx, &packet);

	av_free_packet(&packet);

	if(write_flag)
		av_write_trailer(format_ctx);
	
	avcodec_close(codec_ctx);
	
	if(write_flag)
	{
		avio_closep(&format_ctx->pb);
		avformat_close_input(&format_ctx);
	}
	else
		av_free(codec_ctx);

	sws_freeContext(sws_ctx);
	av_dict_free(&opts);
	av_free(cp_buffer);
//	av_frame_unref(cp_frame);
	av_frame_free(&cp_frame);
	av_free(zoomed_buffer);
//	av_frame_unref(zoomed_frame);
	av_frame_free(&zoomed_frame);
	avformat_free_context(format_ctx);
	return(0);

}

int hue_from_yuv(uint8_t y, uint8_t u, uint8_t v)
{
	int r, g, b;
	int hue;

	if(u>122 && u<134 && v>122 && v<134)//white
		return(361);
	if(y>180 || y<70)
		return(361);
	r = (int) (y + 1.13983 * (v-128));
	g = (int) (y - 0.39465 * (u-128) - 0.58060 * (v-128));
	b = (int) (y + 2.03211 * (u-128));
	if(r==g && g==b)
		return(361);

	if( r>=g && r>=b )
	{
		if(g>=b)
			hue = (int)(60*(g-b)/(r-b));
		else
			hue = (int)(60*(6-(b-g)/(r-g)));
	} 
	else if( g>=r && g>=b )
	{
		if(r>=b)
			hue = (int)(60*(2-(r-b)/(g-b)));
		else
			hue = (int)(60*(2+(b-r)/(g-r)));
	} else
	{
		if(g>=r)
			hue = (int)(60*(4-(g-r)/(b-r)));
		else 
			hue = (int)(60*(4+(r-g)/(b-g)));
	}

	return(hue);
}

int detect_stripe(AVFrame* frame, int x, int y, int w, int h)//return an offset
{
	int i, j;
	int left_anchor = x, right_anchor = x+w;

	for(j=y; j<y+h; j++)
	{
		for(i=x; i<right_anchor; i++)	
		{
			if(*(frame->data[0]+i+j*frame->linesize[0]) > 20 )
			{
				if(i<right_anchor)
				{
					right_anchor = i;
					break;
				}
			}
		}

		for(i=x+w; i>left_anchor; i--)
		{
			if(*(frame->data[0]+i+j*frame->linesize[0]) > 20 )
			{
				if(i>left_anchor)
				{
					left_anchor = i;
					break;
				}
			}		
		}
	}
//	printf("no black, x=%d, w=%d, left=%d, right=%d\n", x, w, left_anchor, right_anchor);

	if( left_anchor==x && right_anchor==x+w )//all black
		return(0);
	else if( left_anchor==x+w && right_anchor==x )//no black or middle black
		return(0);
	else if( left_anchor<x+w && right_anchor>x )//both black
	{
		if( x+w-left_anchor > right_anchor-x &&  x+(left_anchor-w-x)>0 )
			return(left_anchor-x-w+4);
		else if( right_anchor+w < frame->width )
			return(right_anchor-x-4);
		else
			return(0);
	}
	else if( right_anchor == x && x+(left_anchor-w-x)>0 )
		return(left_anchor-w-x+4);//-( w - (left_anchor-x) )
	else if( left_anchor == x+w && right_anchor+w < frame->width )
		return(right_anchor-x-4);
	else
		return(0);
	
		

}


typedef struct {
    /**
     * Number of samples at each PCM value.
     * histogram[0x8000 + i] is the number of samples at value i.
     * The extra element is there for symmetry.
     */
    uint64_t histogram[0x10001];
} VolDetectContext;

#define MAX_DB 91

static inline double logdb(uint64_t v)
{
    double d = v / (double)(0x8000 * 0x8000);
    if (!v)
        return MAX_DB;
    return log(d) * -4.3429448190325182765112891891660508229; /* -10/log(10) */
}


int scale(AVFrame* source,AVFrame** dest,int dest_fmt,int width, int height){
	//uint8_t *dest_buffer;
	//int dest_buffer_size;
	struct SwsContext *sws_ctx;

	(*dest) = av_frame_alloc();
	//dest_buffer_size = avpicture_get_size( dest_fmt, width, height );
	//dest_buffer = (uint8_t *)av_malloc(dest_buffer_size*sizeof(uint8_t));
	//avpicture_fill((AVPicture *)(*dest), dest_buffer, dest_fmt, width, height );
	(*dest)->width = width;
	(*dest)->height = height;
	(*dest)->format = dest_fmt;//frame->format
	av_frame_get_buffer(*dest, 16);
	sws_ctx = sws_getContext(source->width, source->height, source->format, width, height,dest_fmt, SWS_BICUBIC, NULL, NULL, NULL);

	if(sws_scale(sws_ctx, (const uint8_t **)(source->data), source->linesize, 0, source->height, (*dest)->data, (*dest)->linesize)!=(*dest)->height)
	{
		fprintf(stderr, "ERROR in SWS\n");
		sws_freeContext (sws_ctx);
		av_frame_free(dest);
		return 0;
	}
	sws_freeContext(sws_ctx);
	return 1;
}

void initialFrame(AVFrame** dest,int dest_fmt,int width, int height){
	(*dest) = av_frame_alloc();
	(*dest)->width = width;
	(*dest)->height = height;
	(*dest)->format = dest_fmt;//frame->format
	av_frame_get_buffer(*dest, 16);
}

int scale_no_allocated(AVFrame* source,AVFrame* dest,int dest_fmt,int width, int height){
	struct SwsContext *sws_ctx;
	//SWS_BICUBIC
	sws_ctx = sws_getContext(source->width, source->height, source->format, width, height,dest_fmt, SWS_BICUBIC, NULL, NULL, NULL);
    //printf("source width: %d, height: %d\n",source->width,source->height);
    //printf("dst width:  %d, dst height: %d\n ",width,height);
	int check = sws_scale(sws_ctx, (const uint8_t **)(source->data), source->linesize, 0, source->height, dest->data, dest->linesize);
    //printf("check: %d\n", check);
	if(check!=dest->height)
	{
		fprintf(stderr, "ERROR in SWS\n");
		sws_freeContext (sws_ctx);
		return 0;
	}
	sws_freeContext(sws_ctx);
	return 1;
}


int readImg(char* filename,AVFrame** dest){
	AVFormatContext *format_ctx=NULL;
	AVCodecContext *vcodec_ctx=NULL;
	AVCodec *vcodec=NULL;
	AVCodecContext *acodec_ctx=NULL;
	AVCodec *acodec=NULL;
	AVFrame* dec_frame;
	//AVFrame* testFrame;
	int vi=-1;//video stream index
	int ai=-1;//audio stream index
	double fps=0;
	double duration=0;
	int src_width;
	int src_height;
//	int f=0;//count frames
//	double best_intersect=-0.1;
	int stream_finished = 0;
//	int finished = 0;
	init();
	format_ctx = avformat_alloc_context();
	if(open_av(filename, format_ctx, &vcodec_ctx, &acodec_ctx, &vcodec, &acodec, &vi, &ai, &src_width, &src_height, &fps, &duration))
	{	
		fprintf(stderr, "ERROR opening file\n");
		return 0;
	}
//fprintf(stderr,"%d, %d, %d\n",src_width,src_height,fps);
//cout<<src_width<<"\t"<<src_height<<"\t"<<fps<<endl;
	if(src_width<1 || src_height<1 || fps<1 )
	{
		fprintf(stderr, "Invalid width, height, fps or duration\n");
		return 0;
	}

	dec_frame = av_frame_alloc();

	while(!stream_finished)
	{						
		stream_finished = read_av(format_ctx, vcodec_ctx, acodec_ctx, vcodec, acodec, vi, ai, &dec_frame);

		if(dec_frame->pict_type == AV_PICTURE_TYPE_I || dec_frame->pict_type == AV_PICTURE_TYPE_P || dec_frame->pict_type == AV_PICTURE_TYPE_B)//video frame?
		{				
			////to_image(dec_frame,"test.jpg");
			//if(src_width>600||src_height>600){
			//	if(src_width>src_height){
			//		while(src_width>600){
			//			src_width = src_width*2/3;
			//			src_height = src_height*2/3;
			//		}
			//	}
			//	else{
			//		while(src_height>600){
			//			src_width = src_width*2/3;
			//			src_height = src_height*2/3;
			//		}					
			//	}
			//}
			
			
			if(scale(dec_frame,dest,AV_PIX_FMT_YUV420P,src_width, src_height)==0){//AV_PIX_FMT_RGB24
			//if(scale(dec_frame,dest,dec_frame->format,src_width, src_height)==0){//AV_PIX_FMT_RGB24
				avcodec_close(vcodec_ctx);
				if(ai!=-1)//if no audio stream
					avcodec_close(acodec_ctx);
				
				avformat_close_input(&format_ctx);
				av_frame_free(&dec_frame);
				return 0;
			}
			else{
				avcodec_close(vcodec_ctx);
				if(ai!=-1)//if no audio stream
					avcodec_close(acodec_ctx);
				avformat_close_input(&format_ctx);
				av_frame_free(&dec_frame);
				return 1;
			}
		}
	}
	return 0;
}

int fill_iobuffer(void * opaque,uint8_t *buf, int bufsize){
	memcpy(buf,avpacketMemBuff,bufsize*sizeof(char));
	return bufsize;	
}


int readImgFromBuffer(char* buffer, int buffer_size,AVFrame** dest){
	AVFormatContext *format_ctx=NULL;
	AVCodecContext *vcodec_ctx=NULL;
	AVCodec *vcodec=NULL;
	AVCodecContext *acodec_ctx=NULL;
	AVCodec *acodec=NULL;
	AVFrame* dec_frame;
	int vi=-1;//video stream index
	int ai=-1;//audio stream index
	double fps=0;
	double duration=0;
	int src_width;
	int src_height;

	int stream_finished = 0;
	avpacketMemBuff = buffer;

	init();
	
	format_ctx = NULL;
	format_ctx = avformat_alloc_context();
	unsigned char * iobuffer=(unsigned char *)av_malloc(buffer_size);
	AVIOContext *avio =avio_alloc_context(iobuffer, buffer_size,0,NULL,fill_iobuffer,NULL,NULL);
	format_ctx->pb=avio;


	if(open_av("nothing", format_ctx, &vcodec_ctx, &acodec_ctx, &vcodec, &acodec, &vi, &ai, &src_width, &src_height, &fps, &duration))
	{	
		fprintf(stderr, "ERROR opening file\n");
		return 0;
	}
	if(src_width<1 || src_height<1 || fps<1 )
	{
		fprintf(stderr, "Invalid width, height, fps or duration\n");
		return 0;
	}

	dec_frame = av_frame_alloc();



	while(!stream_finished)
	{		
		stream_finished = read_av(format_ctx, vcodec_ctx, acodec_ctx, vcodec, acodec, vi, ai, &dec_frame);
		if(dec_frame->pict_type == AV_PICTURE_TYPE_I || dec_frame->pict_type == AV_PICTURE_TYPE_P || dec_frame->pict_type == AV_PICTURE_TYPE_B)//video frame?
		{				
			if(scale(dec_frame,dest,AV_PIX_FMT_YUV420P,src_width, src_height)==0){//AV_PIX_FMT_RGB24
				if(!(format_ctx->flags & AVFMT_NOFILE)){
					avio_close(format_ctx->pb);
				}
				avcodec_close(vcodec_ctx);
				if(ai!=-1)//if no audio stream
					avcodec_close(acodec_ctx);
				
				avformat_close_input(&format_ctx);
				av_frame_free(&dec_frame);
//				avcodec_close(vcodec);
				return 0;
			}
			else{
				if(!(format_ctx->flags & AVFMT_NOFILE)){
                                        avio_close(format_ctx->pb);
                                }
				avcodec_close(vcodec_ctx);
				if(ai!=-1)//if no audio stream
					avcodec_close(acodec_ctx);
				avformat_close_input(&format_ctx);
				av_frame_free(&dec_frame);
//				avcodec_close(vcodec);
				return 1;
			}
		}
	}
	return 0;
}
