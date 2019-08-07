/*****************************************************************************/
// File: io_examples.cpp
// Author: David Taubman
// Last Revised: 19 August, 2002
/*****************************************************************************/
// Copyright 2001, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"

/*****************************************************************************/
/* EXTERNAL                           main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <threshold> Bi_level\n",argv[0]);
      return -1;
    }

  int T = atoi(argv[3]);

  int err_code=0;
  try {
      int n, width, height, planes;
	  // in is struct of bmp_in defined in io_bmp
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;
	  // width - image columns, height - image rows, planes - number of components (Colour:BGR i.e. 3, Black&White: intensity i.e. 1)
      width = in.cols;  height = in.rows;  planes = in.num_components;
	  // new: similar to malloc, allocate memory; io_byte: 8 bytes, i.e. 0-255 for each variable; dp: data pointer 
      io_byte *dp, *data = new io_byte[width*height*planes];
	  // n=height means read from the bottom row, dp -> data pointer to the allocated memory
	  // for(inital conditions; judgement; operation)
      for (dp=data, n=height; n > 0; n--, dp+=width*planes)
        if ((err_code = bmp_in__get_line(&in,dp)) != 0)
          throw err_code;
      bmp_in__close(&in);

	  // Process the image pixels
	  for (dp = data, n = width*height*planes; n > 0; n--, dp++) {
		  *dp = io_byte((*dp < T)? 0: 255);
	  }

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],width,height,planes)) != 0)
        throw err_code;
      for (dp=data, n=height; n > 0; n--, dp+=width*planes)
        bmp_out__put_line(&out,dp);
      bmp_out__close(&out);

      delete[] data;
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
