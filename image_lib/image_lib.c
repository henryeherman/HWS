/* *********************************************************************************
   EE211A Image Processing Routines (10 October 1999)

   Original Routines written by John Villasenor (villa@icsl.ucla.edu)
   ANSI C/C++ compliant re-write by Matthew Fong (mattfong@icsl.ucla.edu)
 
   Filename:
   image_lib.c

   Description of changes:
   1. Function implementations.
   2. Added DCT image transform that uses the FFT routines.
   3. Fixed error in minimum_value() where zero_flag was used before it was
      initialized. find_zero() initializes this variable.
   4. Pesky data conversion errors when compiled under Visual C++ 6.0 and
      CodeWarrior 5.2 have not been fixed. The problem is not serious, since
      it is a non-explicit double to float conversion. Use gcc if you don't
      want to see it :0). The reason for this is that C++ compilers perform
      very strict type checks. When you multiply two floats together, you get
      a double. The warning occurs when you assign the result to a float.
   5. double to float problems will not be fixed. The FFT2D(), FFT(), DCT2D() 
      and DCT() functions have been changed. The internally allocate double
      rather than float arrays now.
   6. Included <stdlib.h> to prevent compiler errors on Windows machines.
   7. Renames random() to random_number() to prevent conflicts under some
      Linux systems.
   8. Conversion error messages in compliation do not affect program. They
      simply reflect the fact that data has been truncated in the processing
	  of data.

   ****************************************************************************** */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "image_lib.h"

/* #define CHECK_IMAGE_DIM */


/* *************************************************************************** *
                                read_image
 * *************************************************************************** */
void read_image(unsigned char* header, unsigned char** image, char* filename, int n)
{
	FILE *fid;
	int i, j;
	int num_bytes_read = 0;
//	int image_width;
	unsigned char buffer[BUFFER_SIZE];

	if(n > BUFFER_SIZE)
	{
        fprintf(stderr, "[read_byte_array] Image too large. Insufficient buffer size.\n");
        exit(1);	
	}

	/* Open file for binary reading */
	if( (fid = fopen(filename, "rb")) == NULL )
	{
        fprintf(stderr, "[read_byte_array] Cannot open the input file %s.\n", filename);
        exit(1);
	}

	/* Read the header and save it */
	fread(buffer, sizeof(unsigned char), HEADER_LEN, fid);
    for(j = 0; j < HEADER_LEN; j++) 
	{
		header[j] = buffer[j];
	}

#ifdef CHECK_IMAGE_DIM
    /* This code may not always work properly, due to the way ImageMagick writes its 
       .gray files. This code is by default disabled until a better method is found.
       It will be up to the user to determine the correct image size. */
	sscanf( (unsigned char *)(header) + 8, "%7d", &image_width);
	if (image_width != n) 
	{
		fprintf(stderr, "[read_byte_array] Size of image in header (%d) inconsistent with\n", image_width);
		fprintf(stderr, "                  the size n (%d) specified in function call.\n", n);
		exit(1);
	}
#endif
	
	/* Read n bytes at a time */
  	for (i = 0; i < n; i++) 
	{
		num_bytes_read += fread(buffer, sizeof(unsigned char), n, fid);
		
		for (j = 0; j < n; j++) 
		{		
			image[i][j] = buffer[j];
		}
	}

	fclose(fid);
  
	fprintf(stderr, "[read_byte_array] A total of %d bytes were read.\n", num_bytes_read);
}


/* *************************************************************************** *
                                write_image
 * *************************************************************************** */
void write_image(unsigned char* header, unsigned char** image, char* filename, int n)
{
	FILE *fid;
	int i, j;
	int num_bytes_written = 0;
	int total = 0;
	unsigned char buffer[BUFFER_SIZE];

	if(n > BUFFER_SIZE)
	{
        fprintf(stderr, "[write_byte_array] Image too large. Insufficient buffer size.\n");
        exit(1);	
	}

	/* Create file for binary writing */
	if( (fid = fopen(filename, "wb")) == NULL )
	{
        fprintf(stderr, "[write_byte_array] Cannot create the output file %s.\n", filename);
        exit(1);
	}

	/* Write header */
	fwrite(header, sizeof(unsigned char), HEADER_LEN, fid);

	/* Write n bytes at a time */
	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			buffer[j] = image[i][j];
			total += buffer[j];
		}

		num_bytes_written += fwrite(buffer, sizeof(unsigned char), n, fid);
	}
	
	fclose(fid);
  
	fprintf(stderr, "[write_byte_array] A total of %d data bytes were written.\n", num_bytes_written);
	fprintf(stderr, "                   The average value is %d.\n", (total / n) / n);
}


/* *************************************************************************** *
                                write_float_values
 * *************************************************************************** */
void write_float_values(float** image, char* filename, int n)
{
	FILE *fid;
	int i, j;

	/* Create file for binary writing */
	if( (fid = fopen(filename, "wb")) == NULL )
	{
        fprintf(stderr, "[write_float_values] Cannot create the output file %s.\n", filename);
        exit(1);
	}

	/* Write values in row-order */
	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			fprintf(fid, "%f\n", image[i][j]);
		}
	}
	
	fclose(fid);
  
	fprintf(stderr, "[write_float_values] Routine complete.\n");
}


/* *************************************************************************** *
                                generate_header
 * *************************************************************************** */
void generate_header(unsigned char* header, int x_dim, int y_dim)
{
	int i;

	header[0] = '%'; 
	header[1] = 'b'; 
	header[2] = 'i'; 
	header[3] = 't'; 
	header[4] = 'm';
	header[5] = 'a'; 
	header[6] = 'p'; 
	header[7] = '\0'; 

	sprintf( (char *) (header)+8,"%7d", x_dim);
	header[15] = '\0';
	
	sprintf( (char *) (header)+16 ,"%7d", y_dim);
	header[23] = '\0';
	
	/* Number of image planes */
	sprintf( (char *) (header)+24 ,"%7d", 1);
	header[31] = '\0';
	
	/* Number of bits */
	sprintf( (char *) (header)+32 ,"%7d", 8);
	header[39] = '\0';
	
	/* Number of physbits */
	sprintf( (char *) (header)+40 ,"%7d", 8);
	header[47] = '\0';
	sprintf( (char *) (header)+48 ,"%11d", x_dim);
	header[59] = '\0';
	sprintf( (char *) (header)+60 ,"%11d", x_dim * y_dim);
	header[71] = '\0';
	
	/* Mystery variable */
	sprintf( (char *) (header)+72 ,"%11d", 0);
	header[83] = '\0';
	
	/* Aspect ratio */
	sprintf( (char *) (header)+84 ,"%11.6f", 1.0);
	header[95] = '\0';

	for (i = 96; i < 256; i++) header[i] = '\0';
}


/* *************************************************************************** *
                                byte2float
 * *************************************************************************** */
void byte2float(unsigned char** array_in, float** array_out, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			array_out[i][j] = (float)array_in[i][j];
		}
	}
	
	fprintf(stderr, "[byte2float] Routine complete.\n");
}


/* *************************************************************************** *
                                float2complex
 * *************************************************************************** */
void float2complex(float** array_in, Complex** array_out, int n)
{
	int i, j;
   
	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			array_out[i][j].real = array_in[i][j];
			array_out[i][j].imag = 0.0;
		}
	}
  	
	fprintf(stderr, "[float2complex] Routine complete.\n");
}


/* *************************************************************************** *
                                float2byte
 * *************************************************************************** */
void float2byte(float** array_in, unsigned char** array_out, float scale, int n)
{
	int i, j;
	int size = n * n;
	int temp;
	int high_clip = 0;
	int low_clip = 0;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			temp = (int) (array_in[i][j] * scale + 0.5);

			if (temp > 255)
			{
				temp = 255;
				high_clip++;
			}

			if (temp < 0)
			{
				temp = 0;
				low_clip++;
			}
			
			array_out[i][j] = (unsigned char)temp;
		}
	}
  
	fprintf(stderr, "[float2byte] Conversion results:\n");
	fprintf(stderr, "             A total of %d bytes (%5.2f percent) were clipped at 255.\n",
            high_clip, 100.0 * ( (float)high_clip) / ( (float)size) );
	fprintf(stderr, "             A total of %d bytes (%5.2f percent) were clipped at 0.\n",
            low_clip, 100.0 * ( (float)low_clip) / ( (float)size) );
} 


/* *************************************************************************** *
                                array_abs_val
 * *************************************************************************** */
void array_abs_val(Complex** array_in, float** array_out, float* max_val, float* avg_val, int n)
{
	int i, j;
	
	*max_val = -1.0;
	*avg_val = 0.0;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			array_out[i][j] = sqrt(array_in[i][j].real * array_in[i][j].real
			                       + array_in[i][j].imag * array_in[i][j].imag);
			*avg_val += array_out[i][j];
			
			if (array_out[i][j] > *max_val)
			{
				*max_val = array_out[i][j];
			}	
		}
	}
	
	*avg_val /= (float)(n * n);
	
	fprintf(stderr, "[array_abs_val] Results: \n");
	fprintf(stderr, "                The maximum value is %f.\n", *max_val);
	fprintf(stderr, "                The average value is %f.\n", *avg_val);
}


/* *************************************************************************** *
                                FFT2D
   *************************************************************************** *
   Perform a 2D FFT inplace given a complex 2D array. The direction dir is
   used to denote a forward or inverse transform, with 1 denoting forward, and 
   -1 denoting inverse. Only square images with sizes that are powers of 2 can 
   be used!
   *************************************************************************** */
int FFT2D(Complex** image, int n, int direction)
{
	int i, j;
	int m;
	float *real, *imag;

	/* Check that the image is a power of 2 */
	fprintf(stderr, "[FFT2D] Image dimensions are %dx%d.\n", n, n);
	if(!Powerof2(n, &m))
	{
		fprintf(stderr, "[FFT2D] Image size is not a power of 2!\n");
		return(0);
	}

	/* Transform the rows */
	real = (float *)malloc(n * sizeof(float));
	imag = (float *)malloc(n * sizeof(float));
	
	if(real == NULL || imag == NULL)
	{
		fprintf(stderr, "[FFT2D] Memory cannot be allocated.\n");
		return(0);     	
	}

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++) 
		{
			real[j] = image[i][j].real;
			imag[j] = image[i][j].imag;
		}
		
		FFT(direction, m, real, imag);
		
		for(j = 0; j < n; j++) 
		{
			image[i][j].real = real[j];
			image[i][j].imag = imag[j];
		}
	}
   
	free(real);
	free(imag);

	/* Transform the columns */
	real = (float *)malloc(n * sizeof(float));
	imag = (float *)malloc(n * sizeof(float));

	if(real == NULL || imag == NULL)
	{
		fprintf(stderr, "[FFT2D] Memory cannot be allocated.\n");
		return(0);     	
	}

	for(j = 0; j < n; j++) 
	{
		for(i = 0; i < n; i++) 
		{
			real[i] = image[i][j].real;
			imag[i] = image[i][j].imag;
		}
		
		FFT(direction, m, real, imag);
		
		for(i = 0; i < n; i++)
		{
			image[i][j].real = real[i];
			image[i][j].imag = imag[i];
		}
	}
	
	free(real);
	free(imag);

	if (direction == 1)
	{
		fprintf(stderr, "[FFT2D] Forward 2D FFT complete.\n");
	}
	
	if (direction == -1)
	{
		fprintf(stderr, "[FFT2D] Inverse 2D FFT complete.\n");
	}

	return(1);
}


/* *************************************************************************** *
                                FFT
   *************************************************************************** *
   This computes an in-place complex-to-complex FFT x and y are the real and 
   imaginary arrays of 2^m points.
   
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

      Forward formula:
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = Forward transform
              N   /                                n = 0 ... N-1
                  ---
                  k=0

      Inverse formula:
                  N-1
                  ---
                  \          j k 2 pi n / N
      x(k) =       >   X(n) e                    = Inverse transform
                  /                                n = 0 ... N-1
                  ---
                  k=0
 * *************************************************************************** */
int FFT(int dir, int m, float* x, float* y)
{
	long nn,i,i1,j,k,i2,l,l1,l2;
	float c1,c2,tx,ty,t1,t2,u1,u2,z;

	/* Calculate the number of points */
	nn = 1;
	for (i = 0; i < m; i++) nn *= 2;

	/* Do the bit reversal */
	i2 = nn >> 1;
	j = 0;
	for(i = 0; i < nn - 1; i++)
	{
		if (i < j) 
		{
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		
		k = i2;
		
		while (k <= j) 
		{
			j -= k;
			k >>= 1;
		}
		
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	
	for(l = 0; l < m; l++) 
	{
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		
		for (j = 0; j < l1; j++) 
		{
			for (i = j; i < nn; i += l2) 
			{
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1;
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		
		c2 = sqrt((1.0 - c1) / 2.0);
		
		if (dir == 1)
		{
			c2 = -c2;
		}
		
		c1 = sqrt((1.0 + c1) / 2.0);
	}

	return(1);
}


/* *************************************************************************** *
                                DCT2D
   *************************************************************************** */
int DCT2D(float** image, int n, int direction)
{
	float *vector;
	int m;
	int i;
	int j;
	
	/* Check that the image is a power of 2 */
	fprintf(stderr, "[DCT2D] Image dimensions are %dx%d.\n", n, n);
	if(!Powerof2(n, &m))
	{
		fprintf(stderr, "[DCT2D] Image size is not a power of 2!\n");
		return(0);
	}

	/* Transform the rows */
	vector = (float *) malloc(n * sizeof(float));
	
	if(vector == NULL)
	{
		fprintf(stderr, "[DCT2D] Memory cannot be allocated.\n");
		return(0);     	
	}
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			vector[j] = image[i][j];
		}
		
     	if( DCT(direction, m, n, vector) == 0 )
     	{
			fprintf(stderr, "[DCT2D] Memory cannot be allocated in DCT().\n");
			return(0);     	
     	}
     	
     	for(j = 0; j < n; j++)
     	{
     		image[i][j] = vector[j];
     	}
	}
	
	free(vector);
		
	/* Transform the columns */
	vector = (float *) malloc(n * sizeof(float));
	
	if(vector == NULL)
	{
		fprintf(stderr, "[DCT2D] Memory cannot be allocated.\n");
		return(0);     	
	}
	
	for(j = 0; j < n; j++)
	{
		for(i = 0; i < n; i++)
		{
			vector[i] = image[i][j];
		}
		
     	if( DCT(direction, m, n, vector) == 0 )
     	{
			fprintf(stderr, "[DCT2D] Memory cannot be allocated in DCT().\n");
			return(0);     	
     	}
     	
     	for(i = 0; i < n; i++)
     	{
     		image[i][j] = vector[i];
     	}
	}
	
	free(vector);
	
	if (direction == 1)
	{
		fprintf(stderr, "[DCT2D] Forward 2D DCT complete.\n");
	}
	
	if (direction == -1)
	{
		fprintf(stderr, "[DCT2D] Inverse 2D DCT complete.\n");
	}
	
	return(1);
}


/* *************************************************************************** *
                                DCT
   *************************************************************************** */
int DCT(int dir, int m, int n, float* array)
{
	Complex w;
	Complex x;
	Complex *data;
	int i;
	int nhalf;
	float nroot;
	float nfloat;
	float n2float;
	float arg;
	float *real, *imag;

	nfloat = (float) n;
	n2float = nfloat * 2.0;
	nhalf = n / 2;
	nroot= sqrt( ((float) n) / 2.0 );

	/* Take the forward transform */
	if(dir == 1) 
	{
		real = (float *)malloc(n * sizeof(float));
		imag = (float *)malloc(n * sizeof(float));
		data = (Complex *)malloc(n * sizeof(Complex));
		
		if(real == NULL || imag == NULL || data == NULL) return(0);

		for(i = 0; i < nhalf; i++) 
		{
			real[i] = array[2 * i];
			imag[i] = 0.0;
			
			real[i + nhalf] = array[n - 1 - 2 * i];
			imag[i + nhalf] = 0.0;
		}

		FFT(dir, m, real, imag);

		for(i = 0; i < n; i++) 
		{
			arg = -PI * ((float) i) / n2float;
			
			w.real = cos(arg);
			w.imag = sin(arg);
			
			data[i].real = real[i];
			data[i].imag = imag[i];
			
			Complex_multiplication(w, data[i], &x);
			
			array[i] = x.real / nroot;
		}

		array[0] /= sqrt(2.0);

		free(real);
		free(imag);
		free(data);
		
		return(1);
	}

	/* Take the inverse transform */
	if (dir == -1)
	{
		real = (float *)malloc(2 * n * sizeof(float));
		imag = (float *)malloc(2 * n * sizeof(float));
		
		if(real == NULL || imag == NULL) return(0);

		array[0] *= sqrt(2.0);

		for(i = 0; i < n; i++)
		{
			arg = PI * ((float) i) / n2float;
			
			w.real = cos(arg);
			w.imag = sin(arg);
			
			real[i] = w.real * array[i];
			imag[i] = w.imag * array[i];
		}

		real[n] = 0.0;
		imag[n] = 0.0;

		for(i = n + 1; i < 2 * n; i++)
		{
			arg = PI * ((float) i) / n2float;
			
			w.real = -1.0 * cos(arg);
			w.imag = -1.0 * sin(arg);
    
			real[i] = w.real * array[2 * n - i];
			imag[i] = w.imag * array[2 * n - i];
		}

		/* We double the length of n, so we must add one to m */		
		FFT(dir, m + 1, real, imag);

		for(i = 0; i < n; i++)
		{
			array[i] = 0.5 * real[i] / nroot;
		}

		free(real);
		free(imag);
		
		return(1);
	}
}


/* *************************************************************************** *
                                shift_2d_float
 * *************************************************************************** */
void shift_2d_float(float** array, int n)
{
	int row, column;//, index1, index2;
	int nhalf;  

	nhalf = n / 2;

	for(row = 0; row < nhalf; row++)
	{
		for(column = 0; column < nhalf; column++) 
		{
			/* Swap Q2 and Q4 */
			SWAP(array[row][column], array[row + nhalf][column + nhalf]);
			
			/* Swap Q1 and Q3 */
			SWAP(array[row][column + nhalf], array[row + nhalf][column]);
		}
	}
	
	fprintf(stderr, "[shift_2d_float] Routine complete.\n");
}


/* *************************************************************************** *
                                shift_2d_complex
 * *************************************************************************** */
void shift_2d_complex(Complex** array, int n)
{
	int row, column;
	int nhalf;  

	nhalf = n / 2;

	for(row = 0; row < nhalf; row++)
	{
		for(column = 0; column < nhalf; column++) 
		{
			/* Swap Q2 and Q4 */
			SWAP(array[row][column].real, array[row + nhalf][column + nhalf].real);
			SWAP(array[row][column].imag, array[row + nhalf][column + nhalf].imag);
			
			/* Swap Q1 and Q3 */
			SWAP(array[row][column + nhalf].real, array[row + nhalf][column].real);
			SWAP(array[row][column + nhalf].imag, array[row + nhalf][column].imag);
		}
	}
	
	fprintf(stderr, "[shift_2d_complex] Routine complete.\n");
}


/* *************************************************************************** *
                                Powerof2
 * *************************************************************************** */
int Powerof2(int n, int* exponent)
{
	int power = 0;
	int remainder = 0;
	int test = n;
	
	while(test != 1)
	{	
		remainder = test % 2;
		test = test / 2;
		power++;
	}
	
	if(remainder == 0)
	{
		*exponent = power;
		return(1);
	}
	else
	{
		return(0);
	}
}


/* *************************************************************************** *
                                Zonal_LP
 * *************************************************************************** */
void Zonal_LP(Complex** image_complex, float rho_0, int n)
{
	int row, column;
	float filtered_image;
	float rho;
	int nhalf;
	int x, y;
	
	nhalf = n/2;
	
	for(row = 0; row < n; row++)
	{
		for (column = 0; column < n ; column++)
		{
  			x = row - nhalf;
  			y = column - nhalf;
  			rho = sqrt( (float) (x*x + y*y));

   			filtered_image = 0.0;
   			if ( rho <= rho_0)
   			{
   				filtered_image = 1.0;
   			}

    		image_complex[row][column].real *= filtered_image;
    		image_complex[row][column].imag *= filtered_image;
		}
	}
	
	fprintf(stderr, "[Zonal_LP] Routine complete.\n");
}


/* *************************************************************************** *
                                Zonal_HP
 * *************************************************************************** */
void Zonal_HP(Complex** image_complex, float rho_0, int n)
{
	int row, column;
	float filtered_image;
	float rho;
	int nhalf;
	int x, y;
	
	nhalf = n/2;
	
	for(row = 0; row < n; row++)
	{
		for (column = 0; column < n ; column++)
		{
  			x = row - nhalf;
  			y = column - nhalf;
  			rho = sqrt( (float) (x*x + y*y));

   			filtered_image = 0.0;
   			if ( rho >= rho_0)
   			{
   				filtered_image = 1.0;
   			}

    		image_complex[row][column].real *= filtered_image;
    		image_complex[row][column].imag *= filtered_image;
		}
	}
	
	fprintf(stderr, "[Zonal_HP] Routine complete.\n");
}


/* *************************************************************************** *
                                Butterworth_LP
 * *************************************************************************** */
void Butterworth_LP(Complex** image_complex, float rho_0, float order, int n)
{
	int row, column;
	float filtered_image;
	float rho;
	int nhalf;
	int x, y;
	
	nhalf = n/2;
	
	for(row = 0; row < n; row++)
	{
		for (column = 0; column < n ; column++)
		{
  			x = row - nhalf;
  			y = column - nhalf;
  			rho = sqrt( (float) (x*x + y*y));

   			filtered_image = 1.0 / (1.0 + (ROOTTWO - 1.0) * pow(rho / rho_0, 2.0 * order));

    		image_complex[row][column].real *= filtered_image;
    		image_complex[row][column].imag *= filtered_image;
		}
	}
	
	fprintf(stderr, "[Butterworth_LP] Routine complete.\n");
}


/* *************************************************************************** *
                                Butterworth_HP
 * *************************************************************************** */
void Butterworth_HP(Complex** image_complex, float rho_0, float order, int n)
{
	int row, column;
	float filtered_image;
	float rho;
	int nhalf;
	int x, y;
	
	nhalf = n/2;
	
	for(row = 0; row < n; row++)
	{
		for (column = 0; column < n ; column++)
		{
  			x = row - nhalf;
  			y = column - nhalf;
  			rho = sqrt( (float) (x*x + y*y));

   			filtered_image = 1.0 / (1.0 + (ROOTTWO - 1.0) * pow(rho_0 / rho, 2.0 * order));

    		image_complex[row][column].real *= filtered_image;
    		image_complex[row][column].imag *= filtered_image;
		}
	}
	
	fprintf(stderr, "[Butterworth_HP] Routine complete.\n");
}


/* *************************************************************************** *
                                median_filter
 * *************************************************************************** */
void median_filter(float** image_corr, float** med_img, int median_3_or_5, int n)
{
	int rows;
	int cols;
	int i;

	float ** img_tmp;
	/* Dynamically allocate the 2D arrays */
	img_tmp = (float **) malloc(n * sizeof(float *));
	
	for (i = 0; i < n; i++) 
	img_tmp[i] = (float *) malloc(n * sizeof(float));
	
	//float img_tmp[512][512];

	// buffer initialization
	for(rows=0;rows<n;rows++)
		for(cols=0;cols<n;cols++)
			img_tmp[rows][cols]=0;

	if(median_3_or_5 == 3)
	{
		/* Median Filtering on Rows (Length 3) */
		/* Note: We do not perform convolution on the edges. We start at
           (1, 0) and end at (n, n - 1), in order to prevent array overrun. */
		for(rows = 0; rows < n; rows++)
		{
            for(cols = 1; cols < n - 1; cols++)
			{
                img_tmp[rows][cols] = (float)median3( image_corr[rows][cols - 1],
                                                      image_corr[rows][cols],
                                                      image_corr[rows][cols + 1] );
			}
		}

		/* Median Filtering on Columns (Length 3) */
		/* Note: We do not perform convolution on the edges. We start at
           (0, 1) and end at (n - 1, n), in order to prevent array overrun. */
		for(cols = 0; cols < n; cols++)
		{
            for(rows = 1; rows < n - 1; rows++)
			{
                med_img[rows][cols] = (float)median3( img_tmp[rows - 1][cols],
                                                      img_tmp[rows][cols],
                                                      img_tmp[rows + 1][cols] );
			}
		}

	}
	else
	{
        /* Median Filtering on Rows (Length 5) */
        /* Note: We do not perform convolution on the edges. We start at
           (2, 0) and end at (n, n - 2), in order to prevent array overrun */
        for(rows = 0; rows < n; rows++)
		{
            for(cols = 2; cols < n - 2; cols++)
			{
                img_tmp[rows][cols] = (float)median5( image_corr[rows][cols - 2],
                                                      image_corr[rows][cols - 1],
                                                      image_corr[rows][cols],
                                                      image_corr[rows][cols + 1],
                                                      image_corr[rows][cols + 2] );
			}
        }

        /* Median Filtering on Columns (Length 5) */
        /* Note: We do not perform convolution on the edges. We start at
           (0, 2) and end at (n - 2, n), in order to prevent array overrun */
        for(cols = 0; cols < n; cols++)
		{
            for(rows = 2; rows < n - 2; rows++)
			{
                med_img[rows][cols] = (float)median5( img_tmp[rows - 2][cols],
                                                      img_tmp[rows - 1][cols],
                                                      img_tmp[rows][cols],
                                                      img_tmp[rows + 1][cols],
                                                      img_tmp[rows + 2][cols] );
            }
		}
	}

	free(img_tmp);

}


/* *************************************************************************** *
                                min3
 * *************************************************************************** */
float min3(float a, float b, float c)
{
	float x1;
	float x2;

	x1 = MIN(a,b);
	x2 = MIN(b,c);

	return( MIN(x1,x2) );
}


/* *************************************************************************** *
                                median3
 * *************************************************************************** */
float median3(float a, float b, float c)
{
   float x1;

   //x1 = MAX( MIN(a,b), MIN(a,c) );
   //x1 = MAX( x1, MIN(b,c) );

   // 11/24/05 by Hyungjin Kim
   x1 = MAX( (MIN(a,b)), (MIN(a,c)) );
   x1 = MAX( x1, (MIN(b,c)) );

   return x1 ;
}


/* *************************************************************************** *
                                median5
 * *************************************************************************** */
float median5(float a, float b, float c, float d, float e)
{
	float x1;

	x1 = MAX( min3(a,b,c), min3(a,b,d) );
	x1 = MAX( x1, min3(a,b,e) );
	x1 = MAX( x1, min3(a,c,d) );
	x1 = MAX( x1, min3(a,c,e) );
	x1 = MAX( x1, min3(a,d,e) );
	x1 = MAX( x1, min3(a,d,e) );
	x1 = MAX( x1, min3(b,c,d) );
	x1 = MAX( x1, min3(b,c,e) );
	x1 = MAX( x1, min3(b,d,e) );
	x1 = MAX( x1, min3(c,d,e) );

	return x1;
}


/* *************************************************************************** *
                                add_impulse_noise
 * *************************************************************************** *
   Adds impulsive noise. Noise power is uniformly distributed between 0 and 10 
   times ave signal power. percent is the percent of pixels to be replaced by 
   impulsive noise.
 * *************************************************************************** */
void add_impulse_noise(float** image, float percent, int n)
{
	float* noisy_image;
	int corrupt_index;
	int num_corr_points;
	int size = n * n;
	int i;
	int j;
	float avg_signal_power;
	float total = 0.0;

	/* Allocate memory to convert 2D array to 1D */
	noisy_image = (float*) malloc(size * sizeof(float));

	/* Convert 2D array to 1D */
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			noisy_image[(i * n) + j] = image[i][j];
		}
	}
	
	num_corr_points = (int)( (percent / 100.0) * ((float)size) );
	fprintf(stderr, "[add_impulse_noise] Number of corrupt pixels: %d\n", num_corr_points);
   
	/* Get the average signal power */
	for(i = 0; i < size; i++)
	{
		total += noisy_image[i] * noisy_image[i];
	}
	
	avg_signal_power = total / ((float)size);

	/* Add impulsive noise at num_corr_points locations */
	for(i = 0; i < num_corr_points; i++)
	{
		corrupt_index = (int)(uniformly_dist_random_number() * ((float)size));
		noisy_image[corrupt_index] = sqrt(avg_signal_power * 10.0 * uniformly_dist_random_number());
	}
	
	/* Convert 1D array to 2D */
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			image[i][j] = noisy_image[(i * n) + j];
		}
	}
	
	/* Deallocate the memory */
	free(noisy_image);
}


/* *************************************************************************** *
                                add_gaussian_noise
 * *************************************************************************** */
void add_gaussian_noise(float** image, float snr_db, float* mean_noise_power, int n)
{
	float* noisy_image;
	int i;
	int j;
	int size = n * n;
	float avg_signal_power;
	float avg_noise_power;
	float noise_std_dev;
	float total = 0.0;

	/* Allocate memory to convert 2D array to 1D */
	noisy_image = (float*) malloc(size * sizeof(float));

	/* Convert 2D array to 1D */
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			noisy_image[(i * n) + j] = image[i][j];
		}
	}

	/* Get average signal and noise power */
	for(i = 0; i < size; i++)
	{
		total += noisy_image[i] * noisy_image[i];
	}
	
	avg_signal_power = total / ((float)size);
	avg_noise_power = avg_signal_power / ( pow(10.0, snr_db / 10.0) );

	fprintf(stderr, "[add_gaussian_noise] Statistics:");	
	fprintf(stderr, "                     Average signal power: %f\n", avg_signal_power);
	fprintf(stderr, "                     Average noise power: %f\n", avg_noise_power);

	noise_std_dev = sqrt(avg_noise_power);
	
	/* Get the mean of the noise psd. Multiply by fudge factor due to non-ideal pdf's */
	*mean_noise_power = avg_noise_power * (197.0 / 202);

	/* Add corruptive Gaussian noise to the image */
	for(i = 0; i < size; i++)
	{
		noisy_image[i] += gaussian_dist_random_number(noise_std_dev);
	}
	
	/* Convert 1D array to 2D */
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			image[i][j] = noisy_image[(i * n) + j];
		}
	}
	
	/* Deallocate the memory */
	free(noisy_image);
}


/* *************************************************************************** *
                                scale_array
 * *************************************************************************** */
void scale_array(float** image, float scalar, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			image[i][j] *= scalar;
		}
	}
	
	fprintf(stderr, "[scale_array] Routine complete.\n");
}


/* *************************************************************************** *
                                add_scalar
 * *************************************************************************** */
void add_scalar(float** image, float scalar, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			image[i][j] += scalar;
		}
	}
	
	fprintf(stderr, "[add_scalar] Routine complete.\n");
}


/* *************************************************************************** *
                                initialize_float_array
 * *************************************************************************** */
void initialize_float_array(float** image, float value, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			image[i][j] = value;
		}
	}
	
	fprintf(stderr, "[initialize_float_array] Routine complete.\n");
}


/* *************************************************************************** *
                                duplicate_float_array
 * *************************************************************************** */
void duplicate_float_array(float** original, float** copy, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			copy[i][j] = original[i][j];
		}
	}
	
	fprintf(stderr, "[duplicate_float_array] Routine complete.\n");
}


/* *************************************************************************** *
                                array_superposition
 * *************************************************************************** */
void array_superposition(float** x, float** y, float** result, float a, float b, int n)
{
	int i, j;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{
			result[i][j] = (a * x[i][j]) + (b * y[i][j]);
		}
	}
	
	fprintf(stderr, "[array_superposition] Routine complete.\n");
}


/* *************************************************************************** *
                                Complex_addition
 * *************************************************************************** */
void Complex_addition(Complex a, Complex b, Complex* result)
{
	(*result).real = a.real + b.real;
	(*result).imag = a.imag + b.imag;
}


/* *************************************************************************** *
                                Complex_subtration
 * *************************************************************************** */
void Complex_subtration(Complex a, Complex b, Complex* result)
{
	(*result).real = a.real - b.real;
	(*result).imag = a.imag - b.imag;
}


/* *************************************************************************** *
                                Complex_multiplication
 * *************************************************************************** */
void Complex_multiplication(Complex a, Complex b, Complex* result)
{
	(*result).real = a.real * b.real - a.imag * b.imag;
	(*result).imag = a.real * b.imag + a.imag * b.real;
}


/* *************************************************************************** *
                                Complex_division
 * *************************************************************************** */
void Complex_division(Complex a, Complex b, Complex* result)
{
	float u;
	float v;
	float denom;

	if(b.real == 0.0 && b.imag == 0.0)
	{
		fprintf(stderr, "[Complex_division] Divide by zero.\n");
		exit(1);   	
	}

	denom = b.real * b.real + b.imag * b.imag;
	u = b.real / denom;
	v = (-1 * b.imag) / denom;
	
	(*result).real = a.real * u - a.imag * v;
	(*result).real = a.real * v + a.imag * u;
}


/* *************************************************************************** *
                                random_number
 * *************************************************************************** *
   This function implements the minimal random number generator of Park and 
   Miller. It returns a random integer (long) between 0 and 2^31-1 = 2147483647.
   The period of the random sequence is 2147483647.
   Reference: Numerical Recipes in C -- The Art of Scientific Computing, 
              2nd Ed. by William H. Press, William T. Vetterling, 
              Saul A. Teukolsky and Brian P. Flannery.
              Cambridge University Press, Chapter 7, Random Numbers. pp.274-286.
              
   The initial value for random_seed can be anything other than 0.
 * *************************************************************************** */
long random_number(void)
{
	long k;
	static long random_seed = 1;

	k = random_seed / RAN_IQ;
	random_seed = RAN_IA * (random_seed - k * RAN_IQ) - RAN_IR * k;
	
	if (random_seed < 0) random_seed += RAN_IM;
	
	return(random_seed);
}


/* *************************************************************************** *
                                uniformly_dist_random_number
 * *************************************************************************** */
float uniformly_dist_random_number(void)
{
	double x;

	x = (double)random_number();
	
	return((float)(x * RAN_AM));
}


/* *************************************************************************** *
                                gaussian_dist_random_number
 * *************************************************************************** *
   Reference: Numerical Recipes in C -- The Art of Scientific Computing, 
              2nd Ed. by William H. Press, William T. Vetterling, 
              Saul A. Teukolsky and Brian P. Flannery.
              Cambridge University Press, Chapter 7, Random Numbers. pp.287-290.
 * *************************************************************************** */
float gaussian_dist_random_number(float std_dev)
{
	static int iset = 0;
	static float gauss_prev;
	float fac, rsq, v1, v2;

	if(iset == 0)
	{
		do
		{
			v1 = 2.0 * uniformly_dist_random_number() - 1.0;
			v2 = 2.0 * uniformly_dist_random_number() - 1.0;
			
			rsq = (v1 * v1) + (v2 * v2);
		}
		while (rsq >= 1.0 || rsq == 0.0);
		
		fac = sqrt(-2.0 * log(rsq) / rsq) * std_dev;
		gauss_prev = v1 * fac;
		iset = 1;
		
		return(v2 * fac);
	}
	else
	{
		iset = 0;
		return(gauss_prev);
	}
}


/* *************************************************************************** *
                                histogram
 * *************************************************************************** */
void histogram(unsigned char** image, int* hist, int n)
{
	int i;
	int j;
	int size = 256;
   
	for(i = 0; i < size; i++)
	{
		hist[i] = 0;
	}
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			hist[image[i][j]] += 1;    /* Pixel values are between [0, 255] */
		}
	}
}


/* *************************************************************************** *
                                image_statistics
 * *************************************************************************** */
void image_statistics(float** image, float* mean, float* std_dev, int n)
{
	int i;
	int j;
	int size = n * n;
	float total;

	/* Find the mean */   
	total = 0.0;
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			total += image[i][j];
		}
	}
	
	*mean = total / (float)size;
	
	/* Find the standard deviation */   
	total = 0.0;
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			total += (image[i][j] - *mean) * (image[i][j] - *mean);
		}
	}
	
	*std_dev = sqrt(total / (float)size);
}


/* *************************************************************************** *
                                energy_complex
 * *************************************************************************** */
float energy_complex(Complex** image, int n)
{
	int i;
	int j;
	int size = n * n;
	float total;

	/* Find the mean */   
	total = 0.0;
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			total += (image[i][j].real * image[i][j].real) + (image[i][j].imag * image[i][j].imag);
		}
	}
	
	return total;
}


/* *************************************************************************** *
                                nmse
 * *************************************************************************** */
float nmse(float** original, float** reconstructed, int n)
{
	int i;
	int j;
	int size = n * n;
	double total, mean_orig, mean_recon, var_orig, var_recon;

	/* Find the variance of the original image. Assuming zero mean. */   
	mean_orig = 0.0;
	total = 0.0;
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			total += (double)( (original[i][j] - mean_orig) * (original[i][j] - mean_orig) );
		}
	}
	
	var_orig = total / (double)size;

	/* Find the variance of the reconstructed image. Assuming zero mean. */   
	mean_recon = 0.0;
	total = 0.0;
   
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			total += (double)( (original[i][j] - reconstructed[i][j] - mean_recon) * (original[i][j] - reconstructed[i][j] - mean_recon) );
		}
	}
	
	var_recon = total / (double)size;
	
	return((float)(100.0 * var_recon / var_orig));
}


/* *************************************************************************** *
                                maximum_value
 * *************************************************************************** */
float maximum_value(float** image, int positive_or_negative, int n)
{
	int i, j;
	float max_val;

	if(positive_or_negative == 1)
	{
		max_val = -1.0;
	
		for(i = 0; i < n; i++) 
		{
			for(j = 0; j < n; j++) 
			{
				if (image[i][j] > max_val)
				{
					max_val = image[i][j];
				}	
			}
		}
		
		fprintf(stderr, "[maximum_value] Results: \n");
		fprintf(stderr, "                The maximum positive value is %f.\n", max_val);
		
		return max_val;
	}
	else
	{
		max_val = 1.0;
	
		for(i = 0; i < n; i++) 
		{
			for(j = 0; j < n; j++) 
			{
				if (image[i][j] < max_val)
				{
					max_val = image[i][j];
				}	
			}
		}
		
		fprintf(stderr, "[maximum_value] Results: \n");
		fprintf(stderr, "                The maximum negative value is %f.\n", max_val);
		
		return max_val;
	}
}


/* *************************************************************************** *
                                minimum_value
 * *************************************************************************** */
float minimum_value(float** image, int positive_or_negative, int n)
{
	int i, j;
	int zero_flag;
	float min_val;
	float max_val;

	if(positive_or_negative == 1)
	{
		min_val = 4294967296.0;    /* This value is arbitrarily set as 2^32. */
	
		for(i = 0; i < n; i++) 
		{
			for(j = 0; j < n; j++) 
			{				
				if ( (image[i][j] < min_val) && (image[i][j] > 0.0) )
				{
					min_val = image[i][j];
				}	
			}
		}
		
		/* This code allows zero to be a minimum value, if the max_val and min_val are
		   equal and zero is a value in the image. */
		max_val = maximum_value(image, 1, n);
		zero_flag = find_zero(image, n);
		
		if( zero_flag && (max_val == min_val) )
		{
			min_val = 0;
		}
		
		fprintf(stderr, "[minimum_value] Results: \n");
		fprintf(stderr, "                The minimum positive value is %f.\n", min_val);
		
		return min_val;
	}
	else
	{
		min_val = -4294967296.0;    /* This value is arbitrarily set as -2^32. */
	
		for(i = 0; i < n; i++) 
		{
			for(j = 0; j < n; j++) 
			{
				if ( (image[i][j] > min_val) && (image[i][j] < 0.0) )
				{
					min_val = image[i][j];
				}	
			}
		}

		max_val = maximum_value(image, 0, n);
		zero_flag = find_zero(image, n);
		
		if( zero_flag && (max_val == min_val) )
		{
			min_val = 0;
		}
		
		fprintf(stderr, "[minimum_value] Results: \n");
		fprintf(stderr, "                The minimum negative value is %f.\n", min_val);
		
		return min_val;
	}
}


/* *************************************************************************** *
                                find_zero
 * *************************************************************************** */
int find_zero(float** image, int n)
{
	int i, j;
	int zero_flag = 0;

	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
		{				
			if ( image[i][j] == 0.0 )
			{
				zero_flag = 1;
				return zero_flag;
			}	
		}
	}
	
	return zero_flag;
}
