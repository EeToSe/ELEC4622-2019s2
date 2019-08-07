/*****************************************************************************/
// File: image_comps.h
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

/*****************************************************************************/
/* STRUCT                        my_image_comp                               */
/*****************************************************************************/

struct my_image_comp {
	// Data members: (these occupy space in the structure's block of memory)
	int width;
	int height;
	int stride;
	int border; // Extra rows/cols to leave around the boundary
	float *handle; // Points to start of allocated memory buffer
	float *buf; // Points to the first real image sample
	// Function members: (these do not occupy any space in memory)
	my_image_comp()
	{
		width = height = stride = border = 0;  handle = buf = NULL;
	}
	~my_image_comp()
	{
		if (handle != NULL) delete[] handle;
	}
	void init(int height, int width, int border)
	{
		this->width = width;  this->height = height;  this->border = border;
		stride = width + 2 * border;
		if (handle != NULL)
			delete[] handle; // Delete mem allocated by any previous `init' call
		handle = new float[stride*(height + 2 * border)];
		buf = handle + (border*stride) + border;
	}
	void perform_boundary_extension();
	// This function is implemented in "filtering_main.cpp".
};

struct kn_comp {
	// Data members: (these occupy space in the structure's block of memory)
	int width;
	int height;
	int stride;
	int border; // Extra rows/cols to leave around the boundary
	double *handle; // Points to start of allocated memory buffer
	double *buf; // Points to the first real image sample
	// Function members: (these do not occupy any space in memory)
	kn_comp()
	{
		width = height = stride = border = 0;  handle = buf = NULL;
	}
	~kn_comp()
	{
		if (handle != NULL) delete[] handle;
	}
	void init(int height, int width, int border)
	{
		this->width = width;  this->height = height;  this->border = border;
		stride = width + 2 * border;
		if (handle != NULL)
			delete[] handle; // Delete mem allocated by any previous `init' call
		handle = new double[stride*(height + 2 * border)];
		buf = handle + (border*stride) + border;
	}
	// This function is implemented in "filtering_main.cpp".
};

