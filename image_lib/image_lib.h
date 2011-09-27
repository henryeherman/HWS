/* *********************************************************************************
   EE211A Image Processing Routines (11 October 2000)

   Original Routines written by John Villasenor (villa@icsl.ucla.edu)
   ANSI C/C++ compliant re-write by Matthew Fong (mattfong@icsl.ucla.edu)
 
   Filename:
   image_lib.h

   Description of changes:
   1. Function declarations for image_lib.
   2. Added many more constants in.
   3. Problem with the Complex data structure when compliling under Windows.
      Removed the identifier complex after struct to keep with ANSI standard.
   4. Renamed random() to random_number() to prevent conflicts under some
      Linux systems.

   ****************************************************************************** */

#ifndef IMAGE_LIB_H
#define IMAGE_LIB_H


/* Program Constants */
#define HEADER_LEN				256
#define BUFFER_SIZE				1024
#define MAX_NAME				256
#define PMODE					0777
#define ERROR					-1

/* Complex number data structure */
typedef struct 
{
	float real;
	float imag;
} Complex;


/* Macros */
#define SWAP(a, b) { float temp = (a); (a) = (b); (b) = temp; }
#define MAX(x, y) (x > y) ? x : y
#define MIN(x, y) (x < y) ? x : y


/* Function prototypes */
void read_image(unsigned char* header, unsigned char** image, char* filename, int n);
void write_image(unsigned char* header, unsigned char** image, char* filename, int n);
void write_float_values(float** image, char* filename, int n);
void generate_header(unsigned char* header, int x_dim, int y_dim);
void byte2float(unsigned char** array_in, float** array_out, int n);
void float2complex(float** array_in, Complex** array_out, int n);
void float2byte(float** array_in, unsigned char** array_out, float scale, int n);
void array_abs_val(Complex** array_in, float** array_out, float* max_val, float* avg_val, int n);
int FFT2D(Complex** image, int n, int direction);
int FFT(int dir, int m, float* x, float* y);
int DCT2D(float** image, int n, int direction);
int DCT(int dir, int m, int n, float* array);
void shift_2d_float(float** array, int n);
void shift_2d_complex(Complex** array, int n);
int Powerof2(int n, int* exponent);
void Zonal_LP(Complex** image_complex, float rho_0, int n);
void Zonal_HP(Complex** image_complex, float rho_0, int n);
void Butterworth_LP(Complex** image_complex, float rho_0, float order, int n);
void Butterworth_HP(Complex** image_complex, float rho_0, float order, int n);
void median_filter(float** image_corr, float** med_img, int median_3_or_5, int n);
float min3(float a, float b, float c);
float median3(float a, float b, float c);
float median5(float a, float b, float c, float d, float e);
void add_impulse_noise(float** image, float percent, int n);
void add_gaussian_noise(float** image, float snr_db, float* mean_noise_power, int n);
void scale_array(float** image, float scalar, int n);
void add_scalar(float** image, float scalar, int n);
void initialize_float_array(float** image, float value, int n);
void duplicate_float_array(float** original, float** copy, int n); /* How about void** to accept more types ? */
void array_superposition(float** x, float** y, float** result, float a, float b, int n);
void Complex_addition(Complex a, Complex b, Complex* result);
void Complex_subtration(Complex a, Complex b, Complex* result);
void Complex_multiplication(Complex a, Complex b, Complex* result);
void Complex_division(Complex a, Complex b, Complex* result);
long random_number(void);
float uniformly_dist_random_number(void);
float gaussian_dist_random_number(float std_dev);
void histogram(unsigned char** image, int* hist, int n);
void image_statistics(float** image, float* mean, float* std_dev, int n);
float energy_complex(Complex** image, int n);
float nmse(float** original, float** reconstructed, int n);
float maximum_value(float** image, int positive_or_negative, int n);
float minimum_value(float** image, int positive_or_negative, int n);
int find_zero(float** image, int n);


/* Constants needed for random number generation */
#define ROOTTWO					1.4142136
#define RAN_IA					16807
#define RAN_IM					2147483647
#define RAN_IQ					127773
#define RAN_IR					2836
#define RAN_AM					(1.0 / RAN_IM)


/* Numberic constants */
#define PI						3.141592654
#define PIOVER2					1.570796327
#define PIOVER3					1.047197551
#define PIOVER4					0.7853981635
#define TWOPI					6.283185308
#define MINUS_THREE_DB			0.5011872336
#define PLUS_THREE_DB			1.995262315
#define THREE_DB				1.995262315
#define MINUS_FIVE_DB			0.3162277660
#define FIVE_DB					3.162277660
#define PLUS_FIVE_DB			3.162277660
#define MINUS_NINE_DB			0.1258925412
#define NINE_DB					7.943282347
#define PLUS_NINE_DB			7.943282347
#define MINUS_FIFTEEN_DB		0.03162277660
#define FIFTEEN_DB				31.62277660
#define PLUS_FIFTEEN_DB			31.62277660
#define MINUS_TWENTY_DB			0.01         
#define TWENTY_DB				100.0
#define PLUS_TWENTY_DB			100.0
#define MINUS_THIRTY_DB			0.001         
#define THIRTY_DB				1000.0
#define PLUS_THIRTY_DB			1000.0
#define MINUS_FORTY_DB			0.0001         
#define FORTY_DB				10000.0
#define PLUS_FORTY_DB			10000.0


/* Powers of 2 */
#define TWO0					1
#define TWO1					2
#define TWO2					4
#define TWO3					8
#define TWO4					16
#define TWO5					32
#define TWO6					64
#define TWO7					128
#define TWO8					256
#define TWO9					512
#define TWO10					1024
#define TWO11					2048
#define TWO12					4096
#define TWO13					8192
#define TWO14					16384
#define TWO15					32768
#define TWO16					65536
#define TWO17					131072
#define TWO18					262144
#define TWO19					524288
#define TWO20					1048576
#define TWO21					2097152
#define TWO22					4194304
#define TWO23					8388608
#define TWO24					16777216
#define TWO25					33554432
#define TWO26					67108864
#define TWO31					2147483648


/* Physical constants */
#define CLIGHT					2.997924562e8


#endif /* IMAGE_LIB_H */
