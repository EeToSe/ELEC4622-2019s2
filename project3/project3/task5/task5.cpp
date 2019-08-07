/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/

#include "../io_bmp/io_bmp.h"
#include "image_comps.h"
#include <iostream>
#include "motion.h"
#include <math.h>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>


#define vertical 1   
#define horizontal 2
using namespace std;

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*                          Point-symmetry extension			             */
void my_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r * stride + c] = 2 * first_line[c] - first_line[r*stride + c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = 2 * last_line[c] - last_line[-r * stride + c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border * stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = 2 * left_edge[0] - left_edge[c];
			right_edge[c] = 2 * right_edge[0] - right_edge[-c];
		}
}

/*****************************************************************************/
/*           my_image_comp::gaussian_filter   separable convolution         */
/*****************************************************************************/
void gaussian_filter(my_image_comp *in, my_image_comp *out, float sigma, int pattern) {
	constexpr auto PI = 3.141592653589793F;
	int FILTER_EXTENT = ceil(3 * sigma);
	int FILTER_DIM = (2 * FILTER_EXTENT + 1);

	float *filter_buf = new float[FILTER_DIM];
	float *mirror_psf = filter_buf + FILTER_EXTENT;

	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		mirror_psf[t] = exp(-float(t*t) / (2 * sigma*sigma)) / (sigma*sqrt(2 * PI));

	// Check for consistent dimensions
	assert(in->border >= FILTER_EXTENT);
	assert((out->height <= in->height) && (out->width <= in->width));

	// Renormalize to keep the DC gain is 1
	float gain = 0;
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		gain += mirror_psf[t];
	for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
		mirror_psf[t] = mirror_psf[t] / gain;

	if (pattern = vertical) {
		for (int r = 0; r < out->height; r++)
			for (int c = 0; c < out->width; c++) {
				float *ip = in->buf + r * in->stride + c;
				float *op = out->buf + r * out->stride + c;
				float sum = 0.0F;
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					sum += ip[t*in->stride] * mirror_psf[t];
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
				*op = sum;
			}
	}
}

/*****************************************************************************/
/*                     bilinear_interpolation                                */
/*****************************************************************************/
static void
bilinear_interploation(my_image_comp *in, my_image_comp *out, gmvector global_vector) {
	int rr = floor(global_vector.y);
	int cc = floor(global_vector.x);
	float y = global_vector.y - rr;
	float x = global_vector.x - cc;
	float temp;
	for (int r = 0; r < out->height; r++)
		for (int c = 0; c < out->width; c++) {
			float *ip = in->buf + (r - rr) * in->stride + (c - cc);
			float *op = out->buf + r * out->stride + c;
			temp = (1 - x)*(1 - y)*ip[0] + x * (1 - y)*ip[-1] + (1 - x)*y*ip[-in->stride] + x * y*ip[-in->stride - 1];
			*op = int(temp + 0.5);
		}
	cout << "The nearest four integer vectors: ";
	cout << "<" << cc << "," << rr << "> ";
	cout << "<" << (cc + 1) << "," << rr << "> ";
	cout << "<" << cc << "," << (rr + 1) << "> ";
	cout << "<" << (cc + 1) << "," << (rr + 1) << "> " << endl;
	//cout << "y: " << y << "x: " << x;
}

/*****************************************************************************/
/*                              calculate_mse                                */
/*****************************************************************************/
void calculate_mse(my_image_comp *tgt, my_image_comp *output)
{
	float sum = 0;
	for (int r = 0; r < tgt->height; r++) {
		for (int c = 0; c < tgt->width; c++) {
			float *tp = tgt->buf + r * tgt->stride + c;
			float *op = output->buf + r * output->stride + c;
			sum += (tp[0] - op[0]) * (tp[0] - op[0]);
		}
	}
	float mse = sum / (tgt->height*tgt->width);
	std::cout << "The MSE between y'[n] and y[n]: " << mse << endl;
}


/*****************************************************************************/
/* STATIC                         find_motion                                */
/*****************************************************************************/

static mvector
find_motion(my_image_comp *ref, my_image_comp *tgt,
	int start_row, int start_col, int H, int S)
{
	mvector vec, best_vec;
	int sad, best_sad = 256 * (2 * H + 1)*(2 * H + 1);
	// vec.y and vec.x define the search ranges [-8,8] *[-8,8]
	for (vec.y = -S; vec.y <= S; vec.y++)
		for (vec.x = -S; vec.x <= S; vec.x++)
		{
			int ref_row = start_row - vec.y;
			int ref_col = start_col - vec.x;
			if (((ref_row - H) < 0) || ((ref_col - H) < 0) ||
				((ref_row + H) > ref->height) ||
				((ref_col + H) > ref->width))
				continue; // Translated block not containe within reference frame
						  // Continue - not gonna do the following code
			int r, c;
			float *rp = ref->buf + ref_row * ref->stride + ref_col;
			float *tp = tgt->buf + start_row * tgt->stride + start_col;
			for (sad = 0, r = -H; r <= H; r++,
				rp += ref->stride, tp += tgt->stride)
				for (c = -H; c <= H; c++)
				{
					int diff = tp[c] - rp[c];
					sad += (diff < 0) ? (-diff) : diff;  // Absolute it
				}
			if (sad < best_sad)
			{
				best_sad = sad;
				best_vec = vec;
			}
		}
	return best_vec;
}

/*****************************************************************************/
/* STATIC                         motion_comp                                */
/*****************************************************************************/

static void
motion_comp(my_image_comp *ref, my_image_comp *tgt, mvector vec,
	int start_row, int start_col, int block_width, int block_height)
{
	int r, c;
	int ref_row = start_row - vec.y;
	int ref_col = start_col - vec.x;
	float *rp = ref->buf + ref_row * ref->stride + ref_col;
	float *tp = tgt->buf + start_row * tgt->stride + start_col;
	for (r = block_height; r > 0; r--,
		rp += ref->stride, tp += tgt->stride)
		for (c = 0; c < block_width; c++)
			tp[c] = rp[c];
}



/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{
	if (argc != 7)
	{
		fprintf(stderr, "Usage: %s <reference> <target> <out bmp file> <sigma> <block half size: H> <search range: S>\n", argv[0]);
		return -1;
	}
	float sigma = atof(argv[4]);
	int EXTENT = ceil(3 * sigma);


	int H = atoi(argv[5]);
	int S = atoi(argv[6]);
	cout << "***** Command Line Arguments *****\n";
	cout << "Sigma: " << sigma << endl;
	cout << "Block half size: " << H << endl;
	cout << "Search range: " << S << endl;

	int err_code = 0;
	try {
		// Read the input image
		bmp_in in[2];
		if ((err_code = bmp_in__open(&in[0], argv[1])) != 0)
			throw err_code;
		if ((err_code = bmp_in__open(&in[1], argv[2])) != 0)
			throw err_code;

		// Check if the two input frames are the same dimensions
		int width = in[0].cols, height = in[0].rows;
		if ((width != in[1].cols) || (height != in[1].rows))
		{
			fprintf(stderr, "The two input frames have different dimensions.\n");
			return -1;
		}

		// Initialize an array of two class instances
		my_image_comp mono[2];
		mono[0].init(height, width, S);
		mono[1].init(height, width, S);
		int n1 = (width - (H + S + 1));
		int n2 = (height - (H + S + 1));
		cout << "n1: " << n1 << ", n2: " << n2 << endl;

		int n, r, c;
		int num_comps = in[0].num_components; // Monochrome input, i.e. only one colour plane
		io_byte *line = new io_byte[width*num_comps];
		for (n = 0; n < 2; n++)
		{
			for (r = height - 1; r >= 0; r--)
			{ // "r" holds the true row index we are reading, since the image
			  // is stored upside down in the BMP file.
				if ((err_code = bmp_in__get_line(&(in[n]), line)) != 0)
					throw err_code;
				io_byte *src = line; // Points to the first sample of component n
				float *dst = mono[n].buf + r * mono[n].stride;
				for (c = 0; c < width; c++, src += num_comps)
					dst[c] = *src;
			}
			bmp_in__close(&(in[n]));
		}


		// Allocate storage for the filtered output
		my_image_comp *output_comps = new my_image_comp[num_comps];
		my_image_comp *intermedia_comps = new my_image_comp[num_comps];
		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].init(height, width, EXTENT);
			output_comps[n].init(height, width, 0); // Don't need a border for output
		}

		// Process the image, all in floating point (easy)
		// Perform the boundary extension;
		mono[0].perform_boundary_extension();
		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].perform_boundary_extension();
		}

		for (n = 0; n < num_comps; n++)
			gaussian_filter(&(mono[0]), intermedia_comps + n, sigma, vertical);

		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].perform_boundary_extension();
			gaussian_filter(intermedia_comps + n, output_comps + n, sigma, horizontal);
		}
		cout << "Stride: " << output_comps->stride << endl;

		// After we get the gaussian filtered image, try to get the K[n]
		my_image_comp *Kn_comp = new my_image_comp[1];
		Kn_comp[0].init(height, width, 0);
		// Initialize the K[n] to zero 
		for (r = 0; r < height; r++) {
			for (c = 0; c < width; c++) {
				float *op = Kn_comp->buf + r * Kn_comp->stride + c;
				op[0] = 0;
			}
		}

		for (r = H + S; r <= n2; r++) {
			for (c = H + S; c <= n1; c++) {
				float *op = output_comps->buf + r * output_comps->stride + c;
				float *kp = Kn_comp->buf + r * Kn_comp->stride + c;
				float sum_11 = 0, sum_12 = 0, sum_21 = 0, sum_22 = 0;
				for (int Br = -H; Br <= H; Br++)
					for (int Bc = -H; Bc <= H; Bc++) {
						float *dp = op + Br * output_comps->stride + Bc;
						float delta_y = (dp[output_comps->stride] - dp[-output_comps->stride]) / 2;
						float delta_x = (dp[1] - dp[-1]) / 2;
						sum_11 += delta_x * delta_x;
						sum_12 += delta_x * delta_y;
						sum_21 += delta_x * delta_y;
						sum_22 += delta_y * delta_y;
					}
				sum_11 = sum_11;
				sum_12 = sum_12;
				sum_21 = sum_21;
				sum_22 = sum_22;
				float det = sum_11 * sum_22 - sum_12 * sum_21;
				float trace = sum_11 + sum_22;
				if (trace == 0)
					cout << " denominator == 0" << endl;
				kp[0] = det / trace;
			}
		}

		// Find the local maximum
		std::vector<float> Kn_array;
		for (r = H + S; r <= n2; r++)
			for (c = H + S; c <= n1; c++) {
				float *op = Kn_comp->buf + r * Kn_comp->stride + c;
				if ((op[0] > op[1]) && (op[0] > op[-1]) && (op[0] > op[Kn_comp->stride]) && (op[0] > op[-Kn_comp->stride])) {
					Kn_array.push_back(op[0]);
				}
			}
		sort(Kn_array.begin(), Kn_array.end(), greater<float>()); // This the local max vector
		cout << "Kn size: " << Kn_array.size() << endl;

		int Nk;
		cout << "Please enter Nk: " << endl;
		cin >> Nk;
		float tao = Kn_array[Nk];
		cout << "We find tao = " << tao << endl;
		// Free the vector's memory
		Kn_array.clear();
		Kn_array.shrink_to_fit();

		// Find the global vector
		gmvector global_vec;
		int count = 0;
		for (r = H + S; r <= n2; r++)
			for (c = H + S; c <= n1; c++) {
				float *op = Kn_comp->buf + r * Kn_comp->stride + c;
				if ((op[0] > op[1]) && (op[0] > op[-1]) &&
					(op[0] > op[Kn_comp->stride]) && (op[0] > op[-Kn_comp->stride]) && (op[0] > tao)) {
					mvector vec = find_motion(&(mono[0]), &(mono[1]), r, c, H, S);
					global_vec.x += vec.x;
					global_vec.y += vec.y;
					count++;
				}
			}

		global_vec.x = global_vec.x / (count);
		global_vec.y = global_vec.y / (count);
		cout << "Global motion vector: <" << global_vec.x << ", " << global_vec.y << ">\n";

		// Allocate storage for the motion compensated output
		my_image_comp output;
		output.init(height, width, 0); // Don't need a border for output

		// Perform the boundary extension;
		mono[0].perform_boundary_extension();
		// Perform the shifting operation with the bilinear interploation
		bilinear_interploation(&(mono[0]), &output, global_vec);

		calculate_mse(&(mono[1]), &output);
		// Write the motion compensated image out

		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[3], width, height, 1)) != 0)
			throw err_code;
		for (r = height - 1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
		  // written upside down in BMP files.
			io_byte *dst = line; // Points to first sample of component n
			float *src = output.buf + r * output.stride;
			for (int c = 0; c < width; c++, dst++)
				*dst = (io_byte)src[c];
			bmp_out__put_line(&out, line);
		}
		bmp_out__close(&out);
		delete[] line;

		delete[] intermedia_comps;
		delete[] line;
		delete[] output_comps;
		delete[] Kn_comp;
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
