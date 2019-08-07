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
/*              my_aligned_image_comp::perform_boundary_extension            */
/*****************************************************************************/

void my_aligned_image_comp::perform_boundary_extension()
{
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r * stride + c] = first_line[c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border * stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
}


/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
main(int argc, char *argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <half-window size> reduction\n", argv[0]);
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
			input_comps[n].init(height, width, 14); // Leave a border of 4

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

		// Process the image
		int output_width, output_height; // Change the output image size
		output_width = (int)ceil(width * 0.4F);
		output_height = (int)ceil(height * 0.4F);
		// Allocate storage for the filtered output
		my_aligned_image_comp *output_comps =
			new my_aligned_image_comp[num_comps];
		my_aligned_image_comp * intermedia_comps =
			new my_aligned_image_comp[num_comps];
		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].init(output_height, width, 14);  // intermedia image after vertically filtering 
			output_comps[n].init(output_height, output_width, 0); // Don't need a border for output
		}
		for (n = 0; n < num_comps; n++) {
			input_comps[n].perform_boundary_extension();
			intermedia_comps[n].perform_boundary_extension();
		}
		
		// Process the image, all in floating point (easy)
		clock_t start_time = clock();
#if 1	
		cout << "Separable Filtering starts!\n";
		for (n = 0; n < num_comps; n++)
			intermedia_comps[n].filter(input_comps + n, H, vertical);

		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].perform_boundary_extension();
			output_comps[n].filter(intermedia_comps + n, H, horizontal);
		}
		delete[] intermedia_comps;
#elif 0
		cout << "Direct Filtering starts!\n";
		for (n = 0; n < num_comps; n++)
			output_comps[n].filter(input_comps + n, H, direct);
#else 0
		cout << "Vectorized implementation!\n";
		for (n = 0; n < num_comps; n++)
			intermedia_comps[n].vector_filter(input_comps + n, vertical);

		for (n = 0; n < num_comps; n++) {
			intermedia_comps[n].perform_boundary_extension();
			output_comps[n].vector_filter(intermedia_comps + n, horizontal);
		}
		delete[] intermedia_comps;
#endif
		clock_t end_time = clock();
		float elaps = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
		cout << "Filtering done! The computation time: " << elaps << " seconds.\n";


		// Write the image back out again
		bmp_out out;
		if ((err_code = bmp_out__open(&out, argv[2], output_width, output_height, num_comps)) != 0)
			throw err_code;
		io_byte *output_line = new io_byte[output_width*num_comps];  // shrinked line size
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
		delete[] input_comps;
		delete[] output_line;
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
