#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#define vertical 1   
#define horizontal 2
#define direct 3

/*****************************************************************************/
/*       my_aligned_image_comp::filter   separable and 2D convolution        */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int H, int pattern)
{
	constexpr auto PI = 3.141592653589793F;

	int FILTER_EXTENT = H;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);
	float tao = FILTER_EXTENT + 0.5F; // Half length of the hanning window

	// Create the 1D filter PSF as a local array on the stack.
	float *filter_buf0 = new float[FILTER_DIM];  // No shift sinc kernel
	float *filter_buf1 = new float[FILTER_DIM];  // 1/2 shift sinc kernel

	// `mirror_psf' points to the central tap in the filter
	float *mirror_psf0 = filter_buf0 + FILTER_EXTENT;
	float *mirror_psf1 = filter_buf1 + FILTER_EXTENT;

	// Build windowed sinc filter separably
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
	{
		// Special case: the 1-D filter only has one tap
		if (FILTER_EXTENT == 0)
			mirror_psf1[t] = 1;
		else
			mirror_psf1[t] = 0.4F * sinf(0.4F*PI*(t - 0.5F)) / (0.4F*PI*(t - 0.5F)) * 0.5F * (1 + cosf(0.4F*PI*(t - 0.5F) / tao));

		if (t == 0)
			mirror_psf0[t] = 1;
		else
			mirror_psf0[t] = 0.4F * sinf(0.4F*PI*(t)) / (0.4F*PI*t) * 0.5F * (1 + cosf(0.4F*PI*(t) / tao));
	}

	int r, c, rr, cc; // r,c->output , rr,cc->input

	// Create the 2D filter PSF as a local array on the stack.
	if (pattern == direct) {
		int FILTER_TAPS = FILTER_DIM * FILTER_DIM;
		// 0: No shift 2D sinc kernel; 1: Horizontally 1/2 shift 2D sinc kernel; 
		// 2:Vertically 1/2 shift 2D sinc kernel; 3:Vertically and horizontally 1/2 shift 2D sinc kernel
		float **filter_buf = new float *[4];
		float **mirror_psf = new float *[4];

		// Create the 2D filter PSF as a local array on the stack.
		for (int i = 0; i < 4; i++) {
			filter_buf[i] = new float[FILTER_TAPS];
			mirror_psf[i] = filter_buf[i] + (FILTER_DIM*FILTER_EXTENT) + FILTER_EXTENT;
		}

		for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
			for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
				mirror_psf[0][r*FILTER_DIM + c] = mirror_psf0[c] * mirror_psf0[r];
				mirror_psf[1][r*FILTER_DIM + c] = mirror_psf1[c] * mirror_psf0[r];
				mirror_psf[2][r*FILTER_DIM + c] = mirror_psf0[c] * mirror_psf1[r];
				mirror_psf[3][r*FILTER_DIM + c] = mirror_psf1[c] * mirror_psf1[r];
			}

		// Renormalize
		float *gain = new float[4];
		for (int i = 0; i < 4; i++) {
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
					gain[i] += mirror_psf[i][r*FILTER_DIM + c];
		}

		for (int i = 0; i < 4; i++) {
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
					mirror_psf[i][r*FILTER_DIM + c] = mirror_psf[i][r*FILTER_DIM + c] / gain[i];
		}

		// Check for consistent dimensions
		assert(in->border >= FILTER_EXTENT);
		assert((this->height <= in->height) && (this->width <= in->width));

		// Apply four filter kernels to the input image
		for (r = 0, rr = 0; r < height; r += 2, rr += 5)
			for (c = 0, cc = 0; c < width; c += 2, cc += 5)
			{
				float *ip = in->buf + rr * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
						sum += ip[y*in->stride + x] * mirror_psf[0][y*FILTER_DIM + x];
				*op = sum;
			}

		for (r = 0, rr = 0; r < height; r += 2, rr += 5)
			for (c = 1, cc = 2; c < width; c += 2, cc += 5)
			{
				float *ip = in->buf + rr * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
						sum += ip[y*in->stride + x] * mirror_psf[1][y*FILTER_DIM + x];
				*op = sum;
			}

		for (r = 1, rr = 2; r < height; r += 2, rr += 5)
			for (c = 0, cc = 0; c < width; c += 2, cc += 5)
			{
				float *ip = in->buf + rr * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
						sum += ip[y*in->stride + x] * mirror_psf[2][y*FILTER_DIM + x];
				*op = sum;
			}

		for (r = 1, rr = 2; r < height; r += 2, rr += 5)
			for (c = 1, cc = 2; c < width; c += 2, cc += 5)
			{
				float *ip = in->buf + rr * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
						sum += ip[y*in->stride + x] * mirror_psf[3][y*FILTER_DIM + x];
				*op = sum;
			}
	}

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((this->height <= in->height) && (this->width <= in->width));
	// Renormalize to keep the DC gain is 1
	float gain0 = 0, gain1 = 0;
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
	{
		gain0 += mirror_psf0[t];
		gain1 += mirror_psf1[t];
	}

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf0[t] = mirror_psf0[t] / gain0;
		mirror_psf1[t] = mirror_psf1[t] / gain1;
	}

	// Perform 1D vertical convolution, 
	if (pattern == vertical) {
		for (r = 0, rr = 0; r < height; r += 2, rr += 5)
			for (c = 0; c < width; c++)
			{
				float *ip = in->buf + rr * in->stride + c;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t*in->stride] * mirror_psf0[t];
				*op = sum;
			}
		for (r = 1, rr = 2; r < height; r += 2, rr += 5)
			for (c = 0; c < width; c++)
			{
				float *ip = in->buf + rr * in->stride + c;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t*in->stride] * mirror_psf1[t];
				*op = sum;
			}
	}

	// Perform 1D horizaontal convolution  
	if (pattern == horizontal) {
		for (r = 0; r < height; r++) {
			for (c = 0, cc = 0; c < width; c += 2, cc += 5) {
				float *ip = in->buf + r * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t] * mirror_psf0[t];
				*op = sum;
			}

			for (c = 1, cc = 2; c < width; c += 2, cc += 5) {
				float *ip = in->buf + r * in->stride + cc;
				float *op = buf + r * stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t] * mirror_psf1[t];
				*op = sum;
			}
		}
	}
}
