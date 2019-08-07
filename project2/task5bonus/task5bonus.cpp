/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "../task1/io_bmp.h"
#include "../task2/image_comps.h"
#include <iostream>
#include <math.h>
#include <vector>

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
  if (argc != 4) // Must a pair of displacements (i.e. horizontal and vertical)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <radius r> Closing\n",argv[0]);
      return -1;
    }
  
  float raduis = atof(argv[3]);
  std::cout << "radius: " << raduis << "\n";
  int border = raduis + 0.5;
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
		  output_comps[n].init(height, width, border); // Don't need a border for output

	  // Process the image, all in floating point (easy)
	  for (n = 0; n < num_comps; n++)		  
		  input_comps[n].symmetric_extension();

	  // Dynamically assign an array to store all the memory offsets
	  std::vector<int> a_off;
	  std::cout << "Structing set contain these elements" << "\n";
	  for (int y = -border; y <= border; y++) {
		  for (int x = -border; x <= border; x++) {
			  if ((x*x + y * y) < (raduis + 0.5)*(raduis + 0.5)) {
				  std::cout << "<" << y << "," << x << ">  ";
				  // Store in vector of the memory offsets
				  a_off.push_back(y*input_comps[0].stride + x);
			  }
		  }
	  }

	  // Initialize all the pixels of output image to 255
	  for (n = 0; n < num_comps; n++)
		  for (r = 0; r < height; r++)
			  for (c = 0; c < width; c++)
				  output_comps[n].buf[r*output_comps[n].stride + c] = 255;


	  // Closing
	  for (n = 0; n < num_comps; n++) {
		  for (r = 0; r < height; r++) {
			  for (c = 0; c < width; c++) {
				  float *pn = input_comps[n].buf + r * input_comps[n].stride + c;
				  float *op = output_comps[n].buf + r * output_comps[n].stride + c;
				  int val = 0;
				  for (int i = 0; i < a_off.size(); i++) {
					  val = int(pn[a_off[i]]) || val;
				  }
				  if (val == 0) {
					  for (int i = 0; i < a_off.size(); i++)
						  op[a_off[i]] = 0;
				  }
			  }
		  }
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
	// Free the vector's memory
	a_off.clear();
	a_off.shrink_to_fit();
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
