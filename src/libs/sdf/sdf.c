#include <float.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "sdf.h"
void create_sdf(int w, int h, unsigned char *data, struct SDF sdf)
{
    float inv_w = 1.f/(float)w;
    float inv_h = 1.f/(float)h;

    float inv_sdf_w = 1.f/(float)sdf.w;
    float inv_sdf_h = 1.f/(float)sdf.h;
    for(int sdf_y = 0; sdf_y < sdf.h; sdf_y++){
        float sdf_py = (float)sdf_y*inv_sdf_h;
        for(int sdf_x = 0; sdf_x < sdf.w; sdf_x++){
            float sdf_px = (float)sdf_x*inv_sdf_w;
            float min_d2 = FLT_MAX;
            int target_sign;
            {
                int x = (int)(sdf_x*inv_sdf_w*(float)w);
                if(x < 0) x = 0; if(x > w-1) x = w-1;

                int y = (int)(sdf_y*inv_sdf_h*(float)h);
                if(y < 0) y = 0; if(y > h-1) y = h-1;

                target_sign = data[x+y*w];
            }
            for(int y=0;y<h;y++){
                float py = (float)y*inv_h;
                for(int x=0;x<w;x++){
                    if(data[x+y*w] != target_sign){
                        float px = (float)x*inv_w;
                        float dx = px - sdf_px;
                        float dy = py - sdf_py;
                        float d2 = dx*dx + dy*dy;
                        if(d2 < min_d2){
                            min_d2 = d2;
                        }
                    }
                }
            }
            float d = sqrtf(min_d2);
            sdf.data[sdf_x + sdf_y*sdf.w] = d*((float)target_sign*2.f-1.f);
        }
    }

}

unsigned char *get_sdf_file_data(struct SDF sdf, int *len_out)
{
    int len = sizeof(int)*2 + sdf.w*sdf.h*sizeof(float);
    unsigned char *data = calloc(len, sizeof(char));
    unsigned char *d = data;
    memcpy(d,&sdf.w,sizeof(int)); d += sizeof(int);
    memcpy(d,&sdf.h,sizeof(int)); d += sizeof(int);
    memcpy(d,sdf.data,sdf.w*sdf.h*sizeof(float));
    *len_out = len;
    return data;
}

struct SDF sdf_from_file_data(unsigned char *data, int len)
{
    struct SDF ret = {-1,-1,0};
    if(len > sizeof(int)*2){
        int w = *(int*)data; data += sizeof(int);
        int h = *(int*)data; data += sizeof(int);
        if(len >= sizeof(int)*2 + w*h*sizeof(float)){
            ret.w = w;
            ret.h = h;
            ret.data = calloc(w*h, sizeof(float));
            memcpy(ret.data, data, w*h*sizeof(float));
        }
    }
    return ret;
}
