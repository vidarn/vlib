#pragma once
float smootherstep(float edge0, float edge1, float x)
;
float smoothstep(float edge0, float edge1, float x)
;
float clamp(float x, float lowerlimit, float upperlimit)
;
float srgb_to_linear(float srgb)
;
float linear_to_srgb(float linear)
;
void rgb_to_hsv(float* hsv, const float* rgb)
;
void hsv_to_rgb(float* rgb, const float* hsv)
;
void rgb_to_hsl(float* hsl, const float* rgb)
;
void hsl_to_rgb(float* rgb, const float* hsl)
;
void linear_srgb_to_xyz(float* xyz, const float* linear_srgb)
;
void xyz_to_yxy(float* yxy, float* xyz)
;
float yxy_to_cct(float* xyz)
;

