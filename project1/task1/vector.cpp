/*****************************************************************************/
// File: vector_filter.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/
#include <emmintrin.h> // Include SSE2 processor intrinsic functions
#include <stdlib.h>
#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#define vertical 1   
#define horizontal 2

/*****************************************************************************/
/*                     my_aligned_image_comp::vector_filter                  */
/*****************************************************************************/

void my_aligned_image_comp::vector_filter(my_aligned_image_comp *in, int pattern)
{
#define FILTER_EXTENT 0
#define FILTER_TAPS (2*FILTER_EXTENT+1)

	constexpr auto PI = 3.141592653589793F;
	float tao = FILTER_EXTENT + 0.5F;
	// Create the vertical filter PSF as a local array on the stack.
	__m128 filter_buf0[FILTER_TAPS];
	__m128 *mirror_psf0 = filter_buf0 + FILTER_EXTENT;
	__m128 filter_buf1[FILTER_TAPS];
	__m128 *mirror_psf1 = filter_buf1 + FILTER_EXTENT;

	// `mirror_psf' points to the central tap in the filter
	//for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
	//	if (FILTER_EXTENT == 0)
	//		mirror_psf1[t] = _mm_set1_ps(1.0F);
	//	else
	//		mirror_psf1[t] = _mm_set1_ps(0.4F * sinf(PI*(t - 0.5F)) / (PI*(t - 0.5F)) * 0.5F * (1 + cosf(PI*(t - 0.5F) / tao)));

	//	if (t == 0)
	//	mirror_psf0[t] = _mm_set1_ps(1.0F);
	//	else
	//	mirror_psf0[t] = _mm_set1_ps(0.0F);
	//}

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf0[t] = _mm_set1_ps(1.0F / FILTER_TAPS);
		mirror_psf1[t] = _mm_set1_ps(1.0F / FILTER_TAPS);
	}
	__m128 gain0 = _mm_set1_ps(0.0F), gain1 = _mm_set1_ps(0.0F);
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
	{
		gain0 = _mm_add_ps(mirror_psf0[t], gain0);
		gain1 = _mm_add_ps(mirror_psf1[t], gain1);;
	}

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf0[t] = _mm_div_ps(mirror_psf0[t], gain0);
		mirror_psf1[t] = _mm_div_ps(mirror_psf1[t], gain1);
	}
	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((this->height <= in->height) && (this->width <= in->width));
	assert(((stride & 3) == 0) && ((in->stride & 3) == 0));
	int vec_stride_in = in->stride / 4;
	int vec_stride_out = this->stride / 4;
	int vec_width_out = (this->width + 3) / 4; // Big enough to cover the width

	// Do the filtering
	__m128 *out_buf = (__m128 *) buf;
	__m128 *in_buf = (__m128 *)(in->buf);

	if (pattern == vertical) {
		for (int r = 0, rr = 0; r < height; r += 2, rr += 5)
			for (int c = 0; c < vec_width_out; c++)
			{
				__m128 *ip = in_buf + c + rr * vec_stride_in - vec_stride_in * FILTER_EXTENT;
				__m128 *op = out_buf + c + r * vec_stride_out;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++, ip += vec_stride_in)
					sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf0[y], *ip));
				*op = sum;
			}

		for (int r = 1, rr = 2; r < height; r += 2, rr += 5)
			for (int c = 0; c < vec_width_out; c++)
			{
				__m128 *ip = in_buf + c + rr * vec_stride_in - vec_stride_in * FILTER_EXTENT;
				__m128 *op = out_buf + c + r * vec_stride_out;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++, ip += vec_stride_in)
					sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf1[y], *ip));
				*op = sum;
			}
	}

	if (pattern == horizontal) {
		for (int r = 0; r < height; r++) {
			for (int c = 0, cc = 0; c < vec_width_out; c += 2, cc += 5)
			{
				__m128 *ip = in_buf + cc + r * vec_stride_in - FILTER_EXTENT;
				__m128 *op = out_buf + c + r * vec_stride_out;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++, ip+= 1)
					sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf0[y], *ip));
				*op = sum;
			}

			for (int c = 1, cc = 2; c < vec_width_out; c += 2, cc += 5)
			{
				__m128 *ip = in_buf + cc + r * vec_stride_in - FILTER_EXTENT;
				__m128 *op = out_buf + c + r * vec_stride_out;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++, ip+= 1)
					sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf1[y], *ip));
				*op = sum;
			}
		}
	}
}