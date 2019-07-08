#include "common_graphics.h"
#include <math.h>

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
