#include "io_bmp.h"
#include "aligned_image_comps.h"
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;


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
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      first_line[-r*stride+c] = first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r=1; r <= border; r++)
    for (c=0; c < width; c++)
      last_line[r*stride+c] = last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
        left_edge[-c] = left_edge[0];
        right_edge[c] = right_edge[0];
      }
}

void my_aligned_image_comp::compute(my_aligned_image_comp *in1, my_aligned_image_comp *in2)
{
	float sum = 0, sum1 = 0, error = 0;
	for (int r = 0; r < height; r++) {
		for (int c = 0; c < width; c++) {
			float *ip1 = in1->buf + r*in1->stride + c;
			float *ip2 = in2->buf + r*in2->stride + c;
			float *op = this->buf + r*stride + c;
			error = ip1[0] - ip2[0];
			*op = 128 + 0.5F * error;
			sum += error;
			sum1 += pow(error, 2.0F);
		}
	}
	float me = sum / (height*width);
	float mse = sum1 / (height*width);
	float psnr = 10 * log10(255*255/mse);
	cout << "Mean error: " << me << "\nMSE: " << mse << "\nPSNR: " << psnr << "\n"; 
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
      fprintf(stderr,"Usage: %s <in bmp file 0> <in bmp file 1><out bmp file>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in0, in1;
	  if ((err_code = bmp_in__open(&in0, argv[1])) != 0 || (err_code = bmp_in__open(&in1, argv[2])) != 0)
			  throw err_code;

      int width0 = in0.cols, height0 = in0.rows;
	  int width1 = in1.cols, height1 = in1.rows;
      int num_comps = in0.num_components;

      my_aligned_image_comp *input_comps0 =
        new my_aligned_image_comp[num_comps];
	  my_aligned_image_comp *input_comps1 =
		  new my_aligned_image_comp[num_comps];

	  for (int n = 0; n < num_comps; n++) {
		  input_comps0[n].init(height0, width0, 0);
		  input_comps1[n].init(height1, width1, 0);
	  }
      
      int r; // Declare row index
      io_byte *line = new io_byte[width0*num_comps];
      for (r=height0-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in0,line)) != 0)
            throw err_code;
          for (int n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps0[n].buf + r * input_comps0[n].stride;
              for (int c=0; c < width0; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
	  delete[] line;
      bmp_in__close(&in0);

	  io_byte *line1 = new io_byte[width1*num_comps];
	  for (r = height1 - 1; r >= 0; r--)
	  { // "r" holds the true row index we are reading, since the image is
		// stored upside down in the BMP file.
		  if ((err_code = bmp_in__get_line(&in1, line1)) != 0)
			  throw err_code;
		  for (int n = 0; n < num_comps; n++)
		  {
			  io_byte *src = line1 + n; // Points to first sample of component n
			  float *dst = input_comps1[n].buf + r * input_comps1[n].stride;
			  for (int c = 0; c < width1; c++, src += num_comps)
				  dst[c] = (float)*src; 
		  }
	  }
	  delete[] line1;
	  bmp_in__close(&in1);


	  my_aligned_image_comp *output_comps =
		  new my_aligned_image_comp[num_comps];
	  
	  int output_width, output_height;
	  output_width = (width0 < width1) ? width0 : width1;
	  output_height = (height0 < height1) ? height0 : height1;

	  for (int n = 0; n < num_comps; n++) {
		  output_comps[n].init(output_height, output_width, 0); // Don't need a border for output
	  }

	  for (int n = 0; n < num_comps; n++) {
		  cout << "Color plane of " << n + 1 << "\n";
		  output_comps[n].compute(input_comps0 + n, input_comps1 + n);
	  }


      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[3],output_width,output_height,num_comps)) != 0)
        throw err_code;

	  io_byte *output_line = new io_byte[output_width*num_comps];
      for (r=output_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (int n=0; n < num_comps; n++)
            {
              io_byte *dst = output_line+n; // Points to first sample of component n
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
          bmp_out__put_line(&out,output_line);
        }
      bmp_out__close(&out);
      delete[] output_line;
	  delete[] input_comps0;
	  delete[] input_comps1;
      delete[] output_comps;
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
