/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include <iostream>
#include <math.h>


/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::symmetric_extension()
{
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			first_line[-r * stride + c] = first_line[r * stride + c];

	// Now extend downwards
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++)
		for (c = 0; c < width; c++)
			last_line[r*stride + c] = last_line[-r * stride + c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf - border * stride;
	float *right_edge = left_edge + width - 1;
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
		for (c = 1; c <= border; c++)
		{
			left_edge[-c] = left_edge[c];
			right_edge[c] = right_edge[-c];
		}
}

/*****************************************************************************/
/*                              my_image_comp::erosion                       */
/*****************************************************************************/

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if ((argc-1) %2 != 0) // Must a pair of displacements (i.e. horizontal and vertical)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <array of displacements> Erosion\n",argv[0]);
      return -1;
    }
  
  int border = 0; 
  for (int i = 3; i < argc; i++) {
	  int displace = abs(atoi(argv[i]));
	  border = (border> displace)? border:displace;
  }
  std::cout << "border: " << border << "\n";

  int err_code=0;
  try {
	  // Read the input image
	  bmp_in in;
	  if ((err_code = bmp_in__open(&in, argv[1])) != 0)
		  throw err_code;

	  int width = in.cols, height = in.rows;
	  int n, num_comps = in.num_components;
	  my_image_comp *input_comps = new my_image_comp[num_comps];
	  for (n = 0; n < num_comps; n++)
		  input_comps[n].init(height, width, border); // Leave a border of 4

	  std::cout << "stride: " << input_comps[0].stride << "\n";

	  int r,c; // Declare row index
	  io_byte *line = new io_byte[width*num_comps];
	  for (r = height - 1; r >= 0; r--)
	  { // "r" holds the true row index we are reading, since the image is
		// stored upside down in the BMP file.
		  if ((err_code = bmp_in__get_line(&in, line)) != 0)
			  throw err_code;
		  for (n = 0; n < num_comps; n++)
		  {
			  io_byte *src = line + n; // Points to first sample of output image's nth component
			  float *dst = input_comps[n].buf + r * input_comps[n].stride; // Assign the pointer to accordingly input nth component
			  for (int c = 0; c < width; c++, src += num_comps)
				  dst[c] = (float)*src; 
		  }
	  }
	  bmp_in__close(&in);
	 
	  // Allocate storage for the filtered output
	  my_image_comp *output_comps = new my_image_comp[num_comps];
	  for (n = 0; n < num_comps; n++)
		  output_comps[n].init(height, width, 0); // Don't need a border for output

	  // Process the image, all in floating point (easy)
	  for (n = 0; n < num_comps; n++)		  
		  input_comps[n].symmetric_extension();

	  // Create memory offsets
	  int N_A = (argc - 3) / 2;
	  std::cout << N_A << "\n";
	  int *a_off = new int[N_A];
	  for (int i = 3, j = 0; j < N_A; i += 2, j++) 
		a_off[j] = atoi(argv[i])* input_comps[0].stride + atoi(argv[i + 1]);
		  
	  for (int i = 0; i < N_A; i++)
		  std::cout << a_off[i] <<"\n";

	  // Erosion
	  for (n = 0; n < num_comps; n++)
		  for (r = 0; r < height; r++)
			  for (c = 0; c < width; c++)
			  {
				  float *pn = input_comps[n].buf + r*input_comps[n].stride + c;
				  int val = 255;
				  for (int i = 0; i < N_A; i++) {
					  val &= int(pn[a_off[i]]);
				  }
				  output_comps[n].buf[r*output_comps[n].stride + c] = val;
			  }

	  // Write the image back out again
	  bmp_out out;
	  if ((err_code = bmp_out__open(&out, argv[2], width, height, num_comps)) != 0)
		  throw err_code;
	  for (r = height - 1; r >= 0; r--)
	  { // "r" holds the true row index we are writing, since the image is
		// written upside down in BMP files.
		  for (n = 0; n < num_comps; n++)
		  {
			  io_byte *dst = line + n; // Points to first sample of component n
			  float *src = output_comps[n].buf + r * output_comps[n].stride;
			  for (int c = 0; c < width; c++, dst += num_comps) {
				  if (src[c] > 255.05) {
					  src[c] = 255;
				  }
				  else if (src[c] < 0.05) {
					  src[c] = 0;
				  }
				  else
				  {
					  src[c] = int(src[c] + 0.5F);
				  }
				  *dst = (io_byte)src[c];
			  }
		  }
		  bmp_out__put_line(&out, line);
	   }
	bmp_out__close(&out);
	delete[] line;
	delete[] input_comps;
	delete[] output_comps;
	delete[] a_off;
  }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  return 0;
}
