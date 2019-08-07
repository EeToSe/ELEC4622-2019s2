#include<time.h>

#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#define vertical 1
#define horizontal 2
#define num_kernels 5
#define PI 3.141592653589793F
#define SINC(x) (sinf(PI*(x))/(PI*(x)))

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*              my_aligned_image_comp::perform_boundary_extension            */
/*****************************************************************************/

void my_aligned_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r*stride + c] = first_line[c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border*stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}

/*****************************************************************************/
/*                        my_aligned_image_comp::filter                      */
/*****************************************************************************/

void my_aligned_image_comp::filter(my_aligned_image_comp *in, int H, int mode)
{
	int FILTER_EXTENT = H;
	int FILTER_DIM = 2 * FILTER_EXTENT + 1;
	float **filter_buf = new float *[num_kernels];
	float **mirror_psf = new float *[num_kernels];

	// Create the 1D filter PSF as a local array on the stack.
	for (int i = 0; i < num_kernels; i++) {
		filter_buf[i] = new float[FILTER_DIM];
		mirror_psf[i] = filter_buf[i] + FILTER_EXTENT;
	}

	// Build 1D filter kernel
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
		mirror_psf[0][t] = (t == 0)? 1: 0;
		mirror_psf[1][t] = SINC(t - 0.4f)*0.5f*(1 + cosf(PI*(t - 0.4f) / (FILTER_EXTENT + 0.6f)));
		mirror_psf[2][t] = SINC(t + 0.2f)*0.5f*(1 + cosf(PI*(t + 0.2f) / (FILTER_EXTENT + 0.8f))); 
		mirror_psf[3][t] = SINC(t - 0.2f)*0.5f*(1 + cosf(PI*(t - 0.2f) / (FILTER_EXTENT + 0.8f)));
		mirror_psf[4][t] = SINC(t + 0.4f)*0.5f*(1 + cosf(PI*(t + 0.4f) / (FILTER_EXTENT + 0.6f)));
	}

	// Renormalzie
	float *gain = new float[num_kernels];
	for (int i = 0; i < num_kernels; i++) {
		gain[i] = 0;
		for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
			gain[i] += mirror_psf[i][t];
	}

	for (int i = 0; i < num_kernels; i++) 
		for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) 
			mirror_psf[i][t] = mirror_psf[i][t]/gain[i];

	// Perform the convolution
	int r, c, rr, cc; // r,c->output , rr,cc->input
	if (mode == vertical) {
		for (cc = 0; cc < in->width; cc++) {
			for (r = 0, rr = 0; r < height; r += 5, rr += 2) {
				for (int i = 0; i < num_kernels; i++) {
					int nearest_int = int(i*0.4F + 0.5F); // find the nearest interger location
					float *ip = in->buf + ((rr + nearest_int)*in->stride) + cc;
					float *op = this->buf + (r + i)*stride + cc;
					float sum = 0.0F;
					for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
						sum += ip[t*in->stride] * mirror_psf[i][t];
					}
					*op = sum;
				}
			}
		}
	}

	if (mode == horizontal) {
		for (rr = 0; rr < in->height; rr++) {
			for (c = 0, cc = 0; c < width; c += 5, cc += 2) {
				for (int i = 0; i < num_kernels; i++) {
					int nearest_int = int(i*0.4F + 0.5F); // find the nearest interger location
					float *ip = in->buf + rr*in->stride + cc + nearest_int;
					float *op = this->buf + rr*stride + c + i;
					float sum = 0.0F;
					for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++) {
						sum += ip[t] * mirror_psf[i][t];
					}
					*op = sum;
				}
			}
		}
	}
}



/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <half-window size> expansion\n", argv[0]);
		return -1;
	}

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in;
		if ((err_code = bmp_in__open(&in, argv[1])) != 0)
			throw err_code;

		int H = atoi(argv[3]);

		int width = in.cols, height = in.rows;
		int n, num_comps = in.num_components;

		my_aligned_image_comp *input_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++)
			input_comps[n].init(height, width, 8);

		int r; // Declare row index
		io_byte *line = new io_byte[width*num_comps];
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input_comps[n].buf + r * input_comps[n].stride;
				for (int c = 0; c < width; c++, src += num_comps)
					dst[c] = (float)*src; 
			}
		}
		bmp_in__close(&in);

		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermediea_comps =
			new my_aligned_image_comp[num_comps];

		int output_width, output_height;
		output_width = (int)(width*2.5F+0.5F);
		output_height = (int)(height*2.5F+0.5F);

		for (n = 0; n < num_comps; n++) {
			intermediea_comps[n].init(output_height, width, 8);
			output_comps[n].init(output_height, output_width, 0); // Don't need a border for output
		}
		// Process the image, all in floating point (easy)
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
			intermediea_comps[n].perform_boundary_extension();
		}

		for (n = 0; n < num_comps; n++)
			intermediea_comps[n].filter(input_comps + n, H, vertical);
		
		for (n = 0; n < num_comps; n++) {
			intermediea_comps[n].perform_boundary_extension();
			output_comps[n].filter(intermediea_comps + n, H, horizontal);
		}
		delete[] intermediea_comps;

		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], output_width, output_height, num_comps)) != 0)
			throw err_code;

		io_byte *output_line = new io_byte[output_width*num_comps];

		for (r = output_height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			for (n = 0; n < num_comps; n++)
			{
				io_byte *dst = output_line + n; // Points to first sample of component n
				float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < output_width; c++, dst += num_comps) {
					if (src[c] > 255.05) {
						src[c] = 255;
					}
					else if (src[c] < 0.05) {
						src[c] = 0;
					}
					else {
						src[c] = int(src[c] + 0.5F);
					}
					*dst = (io_byte)src[c]; 
				}
			}
			bmp_out__put_line(&out, output_line);
		}
		bmp_out__close(&out);
		delete[] line;
		delete[] output_line;
		delete[] input_comps;
		delete[] output_comps;
	}
	catch (int exc) {
		if (exc == IO_ERR_NO_FILE)
			fprintf(stderr, "Cannot open supplied input or output file.\n");
		else if (exc == IO_ERR_FILE_HEADER)
			fprintf(stderr, "Error encountered while parsing BMP file header.\n");
		else if (exc == IO_ERR_UNSUPPORTED)
			fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
				"simple example supports only 8-bit and 24-bit data.\n");
		else if (exc == IO_ERR_FILE_TRUNC)
			fprintf(stderr, "Input or output file truncated unexpectedly.\n");
		else if (exc == IO_ERR_FILE_NOT_OPEN)
			fprintf(stderr, "Trying to access a file which is not open!(?)\n");
		return -1;
	}
	return 0;
}
