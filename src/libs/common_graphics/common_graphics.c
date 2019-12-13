#include "common_graphics.h"
#include <math.h>

#define min(a,b) ((a < b) ? a : b)
#define max(a,b) ((a > b) ? a : b)

float smootherstep(float edge0, float edge1, float x)
{
	// Scale, and clamp x to 0..1 range
	x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
	// Evaluate polynomial
	return x * x * x * (x * (x * 6 - 15) + 10);
}

float smoothstep(float edge0, float edge1, float x)
{
	// Scale, bias and saturate x to 0..1 range
	x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
	// Evaluate polynomial
	return x * x * (3 - 2 * x);
}

float clamp(float x, float lowerlimit, float upperlimit)
{
	if (x < lowerlimit)
		x = lowerlimit;
	if (x > upperlimit)
		x = upperlimit;
	return x;
}

float srgb_to_linear(float srgb)
{
	if (srgb < 0.04045f) {
		return srgb * (1.f / 12.92f);
	}
	else {
		return powf((srgb + 0.055f) / 1.055f, 2.4f);
	}
}

float linear_to_srgb(float linear)
{
	if (linear < 0.0031308) {
		return linear * 12.92f;
	}
	else {
		return 1.055f*powf(linear, 1.f / 2.4f) - 0.055f;
	}
}

#define MIN_COLOR_DELTA 0.00001f
void rgb_to_hsv(float* hsv, const float* rgb) {
	float min_chroma = min(min(rgb[0], rgb[1]), rgb[2]);
	float max_chroma = max(max(rgb[0], rgb[1]), rgb[2]);
	float delta_chroma = max_chroma - min_chroma;

	hsv[2] = max_chroma;

	if (delta_chroma < MIN_COLOR_DELTA) {
		hsv[0] = 0.f;
		hsv[1] = 0.f;
	}
	else {
		hsv[1] = delta_chroma / max_chroma;

		float delta_r = (((max_chroma - rgb[0]) / 6.f) + (delta_chroma / 2.f)) / delta_chroma;
		float delta_g = (((max_chroma - rgb[1]) / 6.f) + (delta_chroma / 2.f)) / delta_chroma;
		float delta_b = (((max_chroma - rgb[2]) / 6.f) + (delta_chroma / 2.f)) / delta_chroma;

		if (rgb[0] == max_chroma) hsv[0] = delta_b - delta_g;
		else if (rgb[1] == max_chroma) hsv[0] = (1.f / 3.f) + delta_r - delta_b;
		else if (rgb[2] == max_chroma) hsv[0] = (2.f / 3.f) + delta_g - delta_r;

		if (hsv[0] < 0) hsv[0] += 1;
		if (hsv[0] > 1) hsv[0] -= 1;
	}
}
void hsv_to_rgb(float* rgb, const float* hsv) {
	if (hsv[1] == 0.f) {
		rgb[0] = hsv[2];
		rgb[1] = hsv[2];
		rgb[2] = hsv[2];
		return;
	}

	float h = fmodf(hsv[0] * 6.f, 6.f);
	float s = hsv[1];
	float v = hsv[2];

	if (h == 6.f) h = 0.f;
	int h_int = (int)h;
	float var_1 = v * (1 - s);
	float var_2 = v * (1 - s * (h - h_int));
	float var_3 = v * (1 - s * (1 - (h - h_int)));

	float var_r, var_g, var_b;
	if (h_int == 0) { var_r = v; var_g = var_3; var_b = var_1; }
	else if (h_int == 1) { var_r = var_2; var_g = v; var_b = var_1; }
	else if (h_int == 2) { var_r = var_1; var_g = v; var_b = var_3; }
	else if (h_int == 3) { var_r = var_1; var_g = var_2; var_b = v; }
	else if (h_int == 4) { var_r = var_3; var_g = var_1; var_b = v; }
	else { var_r = v; var_g = var_1; var_b = var_2; }

	rgb[0] = var_r;
	rgb[1] = var_g;
	rgb[2] = var_b;
}


static float rotate(float a, float b, float c) {
	if (c < b)
		return a;
	const float delta = c - b;
	while (a < b) a += delta;
	while (a >= c) a -= delta;
	return a;
}
void rgb_to_hsl(float* hsl, const float* rgb) {
	float min_chroma = min(min(rgb[0], rgb[1]), rgb[2]);
	float max_chroma = max(max(rgb[0], rgb[1]), rgb[2]);
	float delta_chroma = max_chroma - min_chroma;

	// Calculate Luminance
	hsl[2] = max((max_chroma + min_chroma)*0.5f, 0.f);

	if (delta_chroma < MIN_COLOR_DELTA) {
		hsl[0] = 0.f;
		hsl[1] = 0.f;
		return;
	}

	// Calculate Hue.
	{
		float h; 
		if (rgb[0] == max_chroma) h = (rgb[1] - rgb[2]) / delta_chroma;
		else if (rgb[1] == max_chroma) h = 2.f + (rgb[2] - rgb[0]) / delta_chroma;
		else if (rgb[2] == max_chroma) h = 4.f + (rgb[0] - rgb[1]) / delta_chroma;

		hsl[0] = rotate(h * 60.f, 0.f, 360.f);

	}

	// Calculate Saturation
	{
		float s;
		// handle the 3 regions defined by [0, 1/2], (1/2, 1], [1, ...]
		if (hsl[2] <= 0.5f)
			s = max((max_chroma - min_chroma) / (max_chroma + min_chroma), 0.0f);
		else if (hsl[2] <= 1.0f)
			s = max((max_chroma - min_chroma) / (2.0f - max_chroma - min_chroma), 0.0f);
		else
			s = max((min_chroma - max_chroma) / (2.0f - max_chroma - min_chroma), 0.0f);

		hsl[1] = s;
	}
}

static float hue_to_rgb(float v1, float v2, float hue) {
	float tmp_h = rotate(hue, 0.f, 1.f);
	if (tmp_h < 1.f/6.f) return v1 + (v2 - v1) * 6.f * tmp_h;
	if (tmp_h < 1.f/2.f) return v2;
	if (tmp_h < 2.f/3.f) return v1 + (v2 - v1) * ((2.f / 3.f) - tmp_h) * 6.f;
	return v1;
}

void hsl_to_rgb(float* rgb, const float* hsl) {
	if (hsl[1] < MIN_COLOR_DELTA) {
		rgb[0] = hsl[2];
		rgb[1] = hsl[2];
		rgb[2] = hsl[2];
		return;
	}
	float h = hsl[0] / 360.f;
	float v1, v2;
	if (hsl[2] <= 0.5f)	v2 = hsl[2] * (1.f + hsl[1]);
	else if (hsl[2] <= 1.0f)	v2 = hsl[2] + hsl[1] - (hsl[2] * hsl[1]);
	else				v2 = hsl[2] - hsl[1] + (hsl[2] * hsl[1]);
	v1 = 2.f * hsl[2] - v2;

	rgb[0] = max(hue_to_rgb(v1, v2, h + 1.f/3.f), 0.f);
	rgb[1] = max(hue_to_rgb(v1, v2, h), 0.f);
	rgb[2] = max(hue_to_rgb(v1, v2, h - 1.f/3.f), 0.f);
}

