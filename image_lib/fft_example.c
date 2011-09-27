/* *********************************************************************************
   EE211A Image Processing Test Program (10 October 1999)
 
   Filename:
   fft_example.c

   Description:
   This program creates a 2D rect and writes it out to a file called rect.gray.
   It is run using the command:

       fft_example output_image_name.gray rect_size image_size output_scale
       
   where filename extension .gray must be included and the image_size must be a
   power of 2.
   
   ****************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image_lib.h"

void main(int argc, char* argv[])
{
	unsigned char** image_byte;
	float** image_float;
	Complex** image_complex;
	unsigned char hdr[HEADER_LEN];
	char output_file[MAX_NAME];
	int rect_size;
	int image_size;
	int output_scale;	
	
	float max_value;
	float avg_value;	
	int i;
	int j;


    /* Absorb command line arguments */
    if ((argc > 5) || (argc < 5))
    {
		fprintf(stderr,
                "Usage: %s [output image name] [rect size] [image size] [output scale]", 
                argv[0]);
        exit(1);
    }
    
    strcpy(output_file, argv[1]);
    rect_size = atoi(argv[2]);
    image_size = atoi(argv[3]);
    output_scale = atoi(argv[4]);    
	
	/* Dynamically allocate the 2D arrays */
	image_byte = (unsigned char **) malloc(image_size * sizeof(unsigned char *));
	image_float = (float **) malloc(image_size * sizeof(float *));
	image_complex = (Complex **) malloc(image_size * sizeof(Complex *));
	
	for (i = 0; i < image_size; i++) 
	{
		image_byte[i] = (unsigned char *) malloc(image_size * sizeof(unsigned char));
		image_float[i] = (float *) malloc(image_size * sizeof(float));
		image_complex[i] = (Complex *) malloc(image_size * sizeof(Complex));
	}

	/* Initialize float array */
	initialize_float_array(image_float, 0.0, image_size);

	/* Define rect function for i, j values between 128 - rect_size/2 and 128 + rect_size/2 */
	for (i = image_size/2 - rect_size/2; i <= image_size/2 + rect_size/2; i++)
	{
		for (j = image_size/2 - rect_size/2; j <= image_size/2 + rect_size/2; j++)
		{
			image_float[i][j] = 255.0;
		}
	}
	
	/* Perform a 2D circular shift */
	shift_2d_float(image_float, image_size);
	
	/* Convert to float to complex */
    float2complex(image_float, image_complex, image_size);

	/* Take the 2D FFT of the complex image */
	if( FFT2D(image_complex, image_size, 1) == 0)
	{
		fprintf(stderr, "Image size is not a power of 2! Exiting.\n");
		exit(1);
	}
	
	/* Take the magnitude of the array */
	array_abs_val(image_complex, image_float, &max_value, &avg_value, image_size);
	
	/* Shift the DC component of the transform to the center of the spectrum */
	shift_2d_float(image_float, image_size);
	
	/* We must scale by 255.0/max_val to get the correct output */
	/* For visibility  purposes, we scale by less */
	float2byte(image_float, image_byte, output_scale/max_value, image_size);

	/* Output to a file */
	generate_header(hdr, image_size, image_size);
	write_image(hdr, image_byte, output_file, image_size);
}
