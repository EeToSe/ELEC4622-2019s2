/*****************************************************************************/
// File: motion_main.cpp
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/

#include "../io_bmp/io_bmp.h"
#include "image_comps.h"
#include "motion.h"
#include <iostream>
using namespace std;

/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/* STATIC                         find_motion                                */
/*****************************************************************************/

static mvector
find_motion(my_image_comp *ref, my_image_comp *tgt,
	int start_row, int start_col, int H, int S)
{
	mvector vec, best_vec;
	int sad, best_sad = 256 * (2*H+1)*(2*H+1);
	// vec.y and vec.x define the search ranges [-8,8] *[-8,8]
	for (vec.y = -S; vec.y <= S; vec.y++)
		for (vec.x = -S; vec.x <= S; vec.x++)
		{
			int ref_row = start_row - vec.y;
			int ref_col = start_col - vec.x;
			if (((ref_row - H) < 0) || ((ref_col - H) < 0) ||
				((ref_row + H) > ref->height) ||
				((ref_col + H) > ref->width)) {
				continue; // Translated block not containe within reference frame
						  // Continue - not gonna do the following code
			}
			int r, c;
			int *rp = ref->buf + ref_row * ref->stride + ref_col;
			int *tp = tgt->buf + start_row * tgt->stride + start_col;
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

/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{
	if (argc != 7)
	{
		fprintf(stderr,
			"Usage: %s <reference frame> <target frame> <bmp MC out> <block half size> <search range> <delta>\n",
			argv[0]);
		return -1;
	}

	int H = atoi(argv[4]);
	int S = atoi(argv[5]);
	int delta = atoi(argv[6]);
	cout << "Block half size: " << H << endl;
	cout << "Search range: " << S << endl;
	cout << "Separation between keypoints: " << delta << endl;


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
		cout << "Width: " << width << ", Height: " << height <<endl;
		mono[0].init(height, width, S); 
		mono[1].init(height, width, S); // Do we need boundary extension?
		int n1 = (width - 2 * (H + S)) / delta; // floor the result
		int n2 = (height - 2 * (H + S)) / delta;
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
				int *dst = mono[n].buf + r * mono[n].stride;
				for (c = 0; c < width; c++, src += num_comps)
					dst[c] = *src;
			}
			bmp_in__close(&(in[n]));
		}

		// Allocate storage for the motion compensated output
		my_image_comp output;
		output.init(height, width, 0); // Don't need a border for output

		// Now perform simple motion estimation and compensation
		gmvector global_vec;
		for (r = 0; r < n2; r ++) {				
			for (c = 0; c < n1; c ++)
			{
				mvector vec = find_motion(&(mono[0]), &(mono[1]),
					(r*delta+H+S), (c*delta+H+S), H, S);
				global_vec.x += vec.x;
				global_vec.y += vec.y;
			}
		}
		global_vec.x = global_vec.x / (n1*n2);
		global_vec.y = global_vec.y / (n1*n2);
		cout << "Global motion vector: <" << global_vec.x << ", " << global_vec.y << ">\n";
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
