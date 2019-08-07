#include <iostream>
#include <cmath>
#include <iomanip>
#include "../io_bmp/io_bmp.h"
#include "image_comps.h"
#include "motion.h"
#include <iostream>
#include <math.h>

using namespace std;
#define vertical 1   
#define horizontal 2

/*****************************************************************************/
/*           my_image_comp::gaussian_filter   separable convolution         */
/*****************************************************************************/
void gaussian_filter(my_image_comp *in, my_image_comp *out, int sigma, int pattern){
	constexpr auto PI = 3.141592653589793F;
	int FILTER_EXTENT = 3*sigma;
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);

	float *filter_buf = new float[FILTER_DIM]; 
	float *mirror_psf = filter_buf + FILTER_EXTENT;
	
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) 
		mirror_psf[t] = exp(-float(t*t)/(2*sigma*sigma))/(sigma*sqrt(2*PI));

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((out->height <= in->height) && (out->width <= in->width));
	
	// Renormalize to keep the DC gain is 1
	float gain = 0;
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		gain += mirror_psf[t];
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		mirror_psf[t] = mirror_psf[t] / gain;
	
	cout << "***** Renormalized Gaussain Filter 1D Kernel ****";
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		cout << mirror_psf[t] << endl;


	if (pattern = vertical) {
		for(int r=0; r < out->height; r++)
			for (int c = 0; c < out->width; c++) {
				float *ip = in->buf + r * in->stride + c;
				float *op = out->buf + r * out->stride + c;
				float sum = 0.0F;
				if ((r < 7) && (c < 7)) {
					cout << "r: "<<r<< " c: "<<c << "value: " << ip[0] <<endl;
				}
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t*in->stride] * mirror_psf[t];
				if (sum < 0) {
					//cout << "****** vertical ******" << endl;
					//cout << sum << "  sum < 0"<<endl;
				}
				*op = sum;
			}
	}
	if (pattern = horizontal) {
		for (int r = 0; r < out->height; r++)
			for (int c = 0; c < out->width; c++) {
				float *ip = in->buf + r * in->stride + c;
				float *op = out->buf + r * out->stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t] * mirror_psf[t];
				if (sum < 0) {
					//cout << "****** horizontal *****" << endl;
					//cout << "sum < 0" << endl;
				}
				*op = sum;
			}
	}
}
