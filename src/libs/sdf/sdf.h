#pragma once
struct SDF
{
    int w;
    int h;
    float *data;
};
void create_sdf(int w, int h, unsigned char *data, struct SDF sdf);
unsigned char *get_sdf_file_data(struct SDF sdf, int *len_out);
struct SDF sdf_from_file_data(unsigned char *data, int len);
