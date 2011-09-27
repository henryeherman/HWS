/* *********************************************************************************
   EE211A Image Processing Test Program (10 October 1999)
 
   Filename:
   gen_rect.c

   Description:
   This program creates a 2D rect and writes it out to a file called rect.gray.
   It is run using the command:

       gen_rect output_image_name.gray rect_size image_size
       
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
	unsigned char hdr[HEADER_LEN];
	char output_file[MAX_NAME];
	int rect_size;
	int image_size;
	
	int i;
	int j;


    /* Absorb command line arguments */
    if ((argc > 4) || (argc < 4))
    {
		fprintf(stderr,
                "Usage: %s [output image name] [rect size] [image size]", 
                argv[0]);
        exit(1);
    }
    
    strcpy(output_file, argv[1]);
    rect_size = atoi(argv[2]);
    image_size = atoi(argv[3]);    
	
	/* Dynamically allocate the 2D arrays */
	image_byte = (unsigned char **) malloc(image_size * sizeof(unsigned char *));
	image_float = (float **) malloc(image_size * sizeof(float *));
	
	for (i = 0; i < image_size; i++) 
	{
		image_byte[i] = (unsigned char *) malloc(image_size * sizeof(unsigned char));
		image_float[i] = (float *) malloc(image_size * sizeof(float));
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
	
	/* Don't forget to shift when you implement the FFT portion of the program! */
	
	/* We must scale by 1.0 to get the correct output */
	float2byte(image_float, image_byte, 1.0, image_size);

	/* Output to a file */
	generate_header(hdr, image_size, image_size);
	write_image(hdr, image_byte, output_file, image_size);	
}
