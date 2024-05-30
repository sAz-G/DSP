/*
 * ----------------------------------------------------------------------------
 * Project Name: fft radix 10 implementation
 * File: fft10.c
 * Description: Impelementation of the radix10 fft function. Originaly implemented to run on a CGRA
 * , although this version deviates a bit from the original version.
 * ----------------------------------------------------------------------------
 * 
 * Author: Sharif Azem
 * GitHub: sAz-G
 * 
 * ----------------------------------------------------------------------------
 * Date Created: 01.06.2023
 * Last Modified: 01.09.2023
 * ----------------------------------------------------------------------------
 * 
 * License: MIT
 * 
 * ----------------------------------------------------------------------------
 * 
 * Change Log:
 * ----------------------------------------------------------------------------
 * 01.06.2023 - implement fft10
 * 01.09.2023 - add simple examples
 * 
 * ----------------------------------------------------------------------------
 * Usage:
 * ----------------------------------------------------------------------------
 * see main function
 * 
 * ----------------------------------------------------------------------------
 * Dependencies:
 * ----------------------------------------------------------------------------
 * 
 * 
 * ----------------------------------------------------------------------------
 * Notes:
 * ----------------------------------------------------------------------------
 * 
 * 
 * ----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_LENGTH 100000000

#ifndef M_PI
#define M_PI           3.142857
#endif

float x[MAX_LENGTH*2]; // array to be transformed 
float X[MAX_LENGTH*2]; // transformed array

void local_fft_radix10(int read_from, int read_jump, int write_to, int lgth, float* x_in, float* X_out);

int main() {
    // Test the local_fft_radix10 function
    int i;
    int length = 100000000; // amount of complex numbers 

    // Initialize the input array x with some values
    for (i = 0; i < length; i++) 
    {    
        x[2*i]   =  i;            // Real part
        x[2*i+1] =  0;            // Imaginary part
    }

    local_fft_radix10(0, 1, 0, length, x, X);

    // Print the output array X
    /*printf("FFT results:\n");
    for (i = 0; i < length*2; i += 2) {
        printf("X[%d] = %g + %gi\n", i/2, X[i], X[i+1]);
    }*/

    /* for (i = 0; i < length*2; i += 2) {
        printf("x[%d] = %.2f + %.2fi\n", i/2, x[i], x[i+1]);
    }*/

    return 0;
}
/*
 * Function: local_fft_radix10
 * ----------------------------------------------------------------------------
 * Description:
 *     Performs a Fast Fourier Transform (FFT) using a radix-10 algorithm.
 * 
 * Parameters:
 *     int read_from:
 *         The starting index from which to read the input data.
 *     int read_jump:
 *         The step size for reading input data.
 *     int write_to:
 *         The starting index where the output data will be written.
 *     int lgth:
 *         The length of the FFT to be performed.
 *     float* x_in:
 *         Pointer to the input array containing the signal to be transformed.
 *     float* X_out:
 *         Pointer to the output array where the transformed signal will be stored.
 * 
 * ----------------------------------------------------------------------------
 * Usage:
 * ----------------------------------------------------------------------------
 * Example call:
 *     int read_from = 0;
 *     int read_jump = 1;
 *     int write_to = 0;
 *     int lgth = 1024;
 *     float x_in[1024];
 *     float X_out[1024];
 *     // Initialize x_in with your data
 *     local_fft_radix10(read_from, read_jump, write_to, lgth, x_in, X_out);
 * 
 * ----------------------------------------------------------------------------
 * Dependencies:
 * ----------------------------------------------------------------------------
 * Standard C Library
 * 
 * ----------------------------------------------------------------------------
 * Notes:
 * ----------------------------------------------------------------------------
 * This implementation assumes that the length of the input (lgth) is a power of 10.
 * Make sure that the input array (x_in) and output array (X_out) are properly allocated.
 * 
 * ----------------------------------------------------------------------------
 */
void local_fft_radix10(int read_from, int read_jump, int write_to, int lgth, float* x_in, float* X_out) {
    if(lgth==1)
    { // copy data to the array to be transformed 
        X_out[2*write_to]   = x_in[2*read_from];
        X_out[2*write_to+1] = x_in[2*read_from+1];
    }
    else
    {
        // even read jump factors  
        local_fft_radix10(   read_from                         ,    10*read_jump,  write_to              ,   lgth/10,   x_in,   X_out  );
        local_fft_radix10(   read_from  +  2*read_jump         ,    10*read_jump,  write_to   + 1*lgth/10,   lgth/10,   x_in,   X_out  );
        local_fft_radix10(   read_from  +  4*read_jump         ,    10*read_jump,  write_to   + 2*lgth/10,   lgth/10,   x_in,   X_out  );
        local_fft_radix10(   read_from  +  6*read_jump         ,    10*read_jump,  write_to   + 3*lgth/10,   lgth/10,   x_in,   X_out  );
        local_fft_radix10(   read_from  +  8*read_jump         ,    10*read_jump,  write_to   + 4*lgth/10,   lgth/10,   x_in,   X_out  );
 
        // odd read jump factors  
        local_fft_radix10(   read_from  +  1*read_jump         ,    10*read_jump,  write_to   + 5*lgth/10,   lgth/10,   x_in,   X_out  );
        local_fft_radix10(   read_from  +  3*read_jump         ,    10*read_jump,  write_to   + 6*lgth/10,   lgth/10,   x_in,   X_out  );         
        local_fft_radix10(   read_from  +  5*read_jump         ,    10*read_jump,  write_to   + 7*lgth/10,   lgth/10,   x_in,   X_out  );         
        local_fft_radix10(   read_from  +  7*read_jump         ,    10*read_jump,  write_to   + 8*lgth/10,   lgth/10,   x_in,   X_out  );         
        local_fft_radix10(   read_from  +  9*read_jump         ,    10*read_jump,  write_to   + 9*lgth/10,   lgth/10,   x_in,   X_out  );         
         
        // for loop to calculate the fft for the given length 
        for(int k=0; k <lgth/10; k++)
        {
            
            int idx_start   = 2*write_to + 2*k; // start reading from here 
            
            // positive frequnecies 
            float even_re_0 = X_out[idx_start];             
            float even_im_0 = X_out[idx_start+1];           

            float even_re_1 = X_out[idx_start + 2*lgth/10];  
            float even_im_1 = X_out[idx_start + 2*lgth/10+1]; 

            float even_re_2 = X_out[idx_start + 4*lgth/10];  
            float even_im_2 = X_out[idx_start + 4*lgth/10+1]; 

            float even_re_3 = X_out[idx_start + 6*lgth/10]; 
            float even_im_3 = X_out[idx_start + 6*lgth/10+1]; 
            
            float even_re_4 = X_out[idx_start + 8*lgth/10]; 
            float even_im_4 = X_out[idx_start + 8*lgth/10+1]; 


            // negative frequencies
            float odd_re_0 = X_out[idx_start+lgth];           
            float odd_im_0 = X_out[idx_start+1+lgth];         

            float odd_re_1 = X_out[idx_start + 2*lgth/10+lgth]; 
            float odd_im_1 = X_out[idx_start + 2*lgth/10+1+lgth];

            float odd_re_2 = X_out[idx_start + 4*lgth/10+lgth]; 
            float odd_im_2 = X_out[idx_start + 4*lgth/10+1+lgth];

            float odd_re_3 = X_out[idx_start + 6*lgth/10+lgth]; 
            float odd_im_3 = X_out[idx_start + 6*lgth/10+1+lgth];

            float odd_re_4 = X_out[idx_start + 8*lgth/10+lgth]; 
            float odd_im_4 = X_out[idx_start + 8*lgth/10+1+lgth];
            
            
            // first frequency set
            float cos_ang00      = cosf(-2.0*M_PI*(1.0/lgth)*k*0)  ;
            float cos_ang01      = cosf(-2.0*M_PI*(1.0/lgth)*k*2.0);
            float cos_ang02      = cosf(-2.0*M_PI*(1.0/lgth)*k*4.0);
            float cos_ang03      = cosf(-2.0*M_PI*(1.0/lgth)*k*6.0);
            float cos_ang04      = cosf(-2.0*M_PI*(1.0/lgth)*k*8.0);

            float sin_ang00      = sinf(-2.0*M_PI*(1.0/lgth)*k*0)  ;
            float sin_ang01      = sinf(-2.0*M_PI*(1.0/lgth)*k*2.0);
            float sin_ang02      = sinf(-2.0*M_PI*(1.0/lgth)*k*4.0);
            float sin_ang03      = sinf(-2.0*M_PI*(1.0/lgth)*k*6.0);
            float sin_ang04      = sinf(-2.0*M_PI*(1.0/lgth)*k*8.0);

            float even_sum_real = even_re_0*cos_ang00    -   even_im_0*sin_ang00 +
                                  even_re_1*cos_ang01    -   even_im_1*sin_ang01 +
                                  even_re_2*cos_ang02    -   even_im_2*sin_ang02 +
                                  even_re_3*cos_ang03    -   even_im_3*sin_ang03 + 
                                  even_re_4*cos_ang04    -   even_im_4*sin_ang04;
                                  
            float even_sum_imag = even_im_0*cos_ang00    +   even_re_0*sin_ang00 +
                                  even_im_1*cos_ang01    +   even_re_1*sin_ang01 +
                                  even_im_2*cos_ang02    +   even_re_2*sin_ang02 +
                                  even_im_3*cos_ang03    +   even_re_3*sin_ang03 +
                                  even_im_4*cos_ang04    +   even_re_4*sin_ang04;

            float odd_sum_real  = odd_re_0*cos_ang00     -   odd_im_0*sin_ang00 +
                                  odd_re_1*cos_ang01     -   odd_im_1*sin_ang01 +
                                  odd_re_2*cos_ang02     -   odd_im_2*sin_ang02 +
                                  odd_re_3*cos_ang03     -   odd_im_3*sin_ang03 +
                                  odd_re_4*cos_ang04     -   odd_im_4*sin_ang04;

            float odd_sum_imag  = odd_im_0*cos_ang00     +   odd_re_0*sin_ang00 +
                                  odd_im_1*cos_ang01     +   odd_re_1*sin_ang01 +
                                  odd_im_2*cos_ang02     +   odd_re_2*sin_ang02 +
                                  odd_im_3*cos_ang03     +   odd_re_3*sin_ang03 +
                                  odd_im_4*cos_ang04     +   odd_re_4*sin_ang04;
           
            float cos_ang_shift0 = cosf(-2.0*M_PI*(1.0/lgth)*k);
            float sin_ang_shift0 = sinf(-2.0*M_PI*(1.0/lgth)*k);


            X_out[idx_start]                 = even_sum_real +  odd_sum_real*cos_ang_shift0  - odd_sum_imag*sin_ang_shift0;
            X_out[idx_start+1]               = even_sum_imag +  odd_sum_imag*cos_ang_shift0  + odd_sum_real*sin_ang_shift0;   
            X_out[idx_start + lgth]          = even_sum_real - (odd_sum_real*cos_ang_shift0 - odd_sum_imag*sin_ang_shift0);
            X_out[idx_start + 1 + lgth]      = even_sum_imag - (odd_sum_imag*cos_ang_shift0 + odd_sum_real*sin_ang_shift0);



            // second frequency set
            float cos_ang10      = cosf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*0);
            float cos_ang11      = cosf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*2);
            float cos_ang12      = cosf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*4);
            float cos_ang13      = cosf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*6);
            float cos_ang14      = cosf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*8);

            float sin_ang10      = sinf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*0);
            float sin_ang11      = sinf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*2);
            float sin_ang12      = sinf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*4);
            float sin_ang13      = sinf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*6);
            float sin_ang14      = sinf(-2.0*M_PI*(1.0/lgth)*(k+lgth/10)*8);

            even_sum_real =     even_re_0*cos_ang10    -   even_im_0*sin_ang10 +
                                even_re_1*cos_ang11    -   even_im_1*sin_ang11 +
                                even_re_2*cos_ang12    -   even_im_2*sin_ang12 +
                                even_re_3*cos_ang13    -   even_im_3*sin_ang13 + 
                                even_re_4*cos_ang14    -   even_im_4*sin_ang14;

            even_sum_imag =     even_im_0*cos_ang10    +   even_re_0*sin_ang10 +
                                even_im_1*cos_ang11    +   even_re_1*sin_ang11 +
                                even_im_2*cos_ang12    +   even_re_2*sin_ang12 +
                                even_im_3*cos_ang13    +   even_re_3*sin_ang13 +
                                even_im_4*cos_ang14    +   even_re_4*sin_ang14;

            odd_sum_real  =     odd_re_0*cos_ang10     -   odd_im_0*sin_ang10 +
                                odd_re_1*cos_ang11     -   odd_im_1*sin_ang11 +
                                odd_re_2*cos_ang12     -   odd_im_2*sin_ang12 +
                                odd_re_3*cos_ang13     -   odd_im_3*sin_ang13 +
                                odd_re_4*cos_ang14     -   odd_im_4*sin_ang14;

            odd_sum_imag  =     odd_im_0*cos_ang10     +   odd_re_0*sin_ang10 +
                                odd_im_1*cos_ang11     +   odd_re_1*sin_ang11 +
                                odd_im_2*cos_ang12     +   odd_re_2*sin_ang12 +
                                odd_im_3*cos_ang13     +   odd_re_3*sin_ang13 +
                                odd_im_4*cos_ang14     +   odd_re_4*sin_ang14;

            float cos_ang_shift1 = cosf(-2.0*M_PI*(1.0/lgth)*(k+1*lgth/10));
            float sin_ang_shift1 = sinf(-2.0*M_PI*(1.0/lgth)*(k+1*lgth/10));


            X_out[idx_start + 2*lgth/10]          = even_sum_real +  odd_sum_real*cos_ang_shift1  - odd_sum_imag*sin_ang_shift1;        
            X_out[idx_start + 2*lgth/10+1]        = even_sum_imag +  odd_sum_imag*cos_ang_shift1  + odd_sum_real*sin_ang_shift1;       
            X_out[idx_start + 2*lgth/10 + lgth]   = even_sum_real - (odd_sum_real*cos_ang_shift1  - odd_sum_imag*sin_ang_shift1);
            X_out[idx_start + 2*lgth/10+1 + lgth] = even_sum_imag - (odd_sum_imag*cos_ang_shift1  + odd_sum_real*sin_ang_shift1);




            // third frequency set
            float cos_ang20      = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*0);
            float cos_ang21      = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*2);
            float cos_ang22      = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*4);
            float cos_ang23      = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*6);
            float cos_ang24      = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*8);

            float sin_ang20      = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*0);
            float sin_ang21      = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*2);
            float sin_ang22      = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*4);
            float sin_ang23      = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*6);
            float sin_ang24      = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10)*8);

            even_sum_real =     even_re_0*cos_ang20    -   even_im_0*sin_ang20 +
                                even_re_1*cos_ang21    -   even_im_1*sin_ang21 +
                                even_re_2*cos_ang22    -   even_im_2*sin_ang22 +
                                even_re_3*cos_ang23    -   even_im_3*sin_ang23 + 
                                even_re_4*cos_ang24    -   even_im_4*sin_ang24;

            even_sum_imag =     even_im_0*cos_ang20    +   even_re_0*sin_ang20 +
                                even_im_1*cos_ang21    +   even_re_1*sin_ang21 +
                                even_im_2*cos_ang22    +   even_re_2*sin_ang22 +
                                even_im_3*cos_ang23    +   even_re_3*sin_ang23 +
                                even_im_4*cos_ang24    +   even_re_4*sin_ang24;

            odd_sum_real  =     odd_re_0*cos_ang20     -   odd_im_0*sin_ang20 +
                                odd_re_1*cos_ang21     -   odd_im_1*sin_ang21 +
                                odd_re_2*cos_ang22     -   odd_im_2*sin_ang22 +
                                odd_re_3*cos_ang23     -   odd_im_3*sin_ang23 +
                                odd_re_4*cos_ang24     -   odd_im_4*sin_ang24;

            odd_sum_imag  =     odd_im_0*cos_ang20     +   odd_re_0*sin_ang20 +
                                odd_im_1*cos_ang21     +   odd_re_1*sin_ang21 +
                                odd_im_2*cos_ang22     +   odd_re_2*sin_ang22 +
                                odd_im_3*cos_ang23     +   odd_re_3*sin_ang23 +
                                odd_im_4*cos_ang24     +   odd_re_4*sin_ang24;     

            float cos_ang_shift2 = cosf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10));
            float sin_ang_shift2 = sinf(-2.0*M_PI*(1.0/lgth)*(k+2*lgth/10));


            X_out[idx_start + 4*lgth/10]          = even_sum_real +  odd_sum_real*cos_ang_shift2  - odd_sum_imag*sin_ang_shift2;
            X_out[idx_start + 4*lgth/10+1]        = even_sum_imag +  odd_sum_imag*cos_ang_shift2  + odd_sum_real*sin_ang_shift2;
            X_out[idx_start + 4*lgth/10 + lgth]   = even_sum_real - (odd_sum_real*cos_ang_shift2  - odd_sum_imag*sin_ang_shift2);
            X_out[idx_start + 4*lgth/10+1 + lgth] = even_sum_imag - (odd_sum_imag*cos_ang_shift2  + odd_sum_real*sin_ang_shift2);




            // fourh frequency set 
            float cos_ang30      = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*0);
            float cos_ang31      = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*2);
            float cos_ang32      = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*4);
            float cos_ang33      = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*6);
            float cos_ang34      = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*8);

            float sin_ang30      = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*0);
            float sin_ang31      = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*2);
            float sin_ang32      = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*4);
            float sin_ang33      = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*6);
            float sin_ang34      = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10)*8);

            even_sum_real =     even_re_0*cos_ang30    -   even_im_0*sin_ang30 +
                                even_re_1*cos_ang31    -   even_im_1*sin_ang31 +
                                even_re_2*cos_ang32    -   even_im_2*sin_ang32 +
                                even_re_3*cos_ang33    -   even_im_3*sin_ang33 + 
                                even_re_4*cos_ang34    -   even_im_4*sin_ang34;

            even_sum_imag =     even_im_0*cos_ang30    +   even_re_0*sin_ang30 +
                                even_im_1*cos_ang31    +   even_re_1*sin_ang31 +
                                even_im_2*cos_ang32    +   even_re_2*sin_ang32 +
                                even_im_3*cos_ang33    +   even_re_3*sin_ang33 +
                                even_im_4*cos_ang34    +   even_re_4*sin_ang34;

            odd_sum_real  =     odd_re_0*cos_ang30     -   odd_im_0*sin_ang30 +
                                odd_re_1*cos_ang31     -   odd_im_1*sin_ang31 +
                                odd_re_2*cos_ang32     -   odd_im_2*sin_ang32 +
                                odd_re_3*cos_ang33     -   odd_im_3*sin_ang33 +
                                odd_re_4*cos_ang34     -   odd_im_4*sin_ang34;

            odd_sum_imag  =     odd_im_0*cos_ang30     +   odd_re_0*sin_ang30 +
                                odd_im_1*cos_ang31     +   odd_re_1*sin_ang31 +
                                odd_im_2*cos_ang32     +   odd_re_2*sin_ang32 +
                                odd_im_3*cos_ang33     +   odd_re_3*sin_ang33 +
                                odd_im_4*cos_ang34     +   odd_re_4*sin_ang34;     


            float cos_ang_shift3 = cosf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10));
            float sin_ang_shift3 = sinf(-2.0*M_PI*(1.0/lgth)*(k+3*lgth/10));


            X_out[idx_start + 6*lgth/10]            = even_sum_real + odd_sum_real*cos_ang_shift3  - odd_sum_imag*sin_ang_shift3;
            X_out[idx_start + 6*lgth/10+1]          = even_sum_imag + odd_sum_imag*cos_ang_shift3  + odd_sum_real*sin_ang_shift3;
            X_out[idx_start + 6*lgth/10 + lgth]     = even_sum_real - (odd_sum_real*cos_ang_shift3 - odd_sum_imag*sin_ang_shift3);
            X_out[idx_start + 6*lgth/10+1 + lgth]   = even_sum_imag - (odd_sum_imag*cos_ang_shift3 + odd_sum_real*sin_ang_shift3);




            // fifth frequency set
            float cos_ang40      = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*0);
            float cos_ang41      = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*2);
            float cos_ang42      = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*4);
            float cos_ang43      = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*6);
            float cos_ang44      = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*8);

            float sin_ang40      = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*0);
            float sin_ang41      = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*2);
            float sin_ang42      = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*4);
            float sin_ang43      = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*6);
            float sin_ang44      = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10)*8);

            even_sum_real =     even_re_0*cos_ang40    -   even_im_0*sin_ang40 +
                                even_re_1*cos_ang41    -   even_im_1*sin_ang41 +
                                even_re_2*cos_ang42    -   even_im_2*sin_ang42 +
                                even_re_3*cos_ang43    -   even_im_3*sin_ang43 + 
                                even_re_4*cos_ang44    -   even_im_4*sin_ang44;

            even_sum_imag =     even_im_0*cos_ang40    +   even_re_0*sin_ang40 +
                                even_im_1*cos_ang41    +   even_re_1*sin_ang41 +
                                even_im_2*cos_ang42    +   even_re_2*sin_ang42 +
                                even_im_3*cos_ang43    +   even_re_3*sin_ang43 +
                                even_im_4*cos_ang44    +   even_re_4*sin_ang44;

            odd_sum_real  =     odd_re_0*cos_ang40     -   odd_im_0*sin_ang40 +
                                odd_re_1*cos_ang41     -   odd_im_1*sin_ang41 +
                                odd_re_2*cos_ang42     -   odd_im_2*sin_ang42 +
                                odd_re_3*cos_ang43     -   odd_im_3*sin_ang43 +
                                odd_re_4*cos_ang44     -   odd_im_4*sin_ang44;

            odd_sum_imag  =     odd_im_0*cos_ang40     +   odd_re_0*sin_ang40 +
                                odd_im_1*cos_ang41     +   odd_re_1*sin_ang41 +
                                odd_im_2*cos_ang42     +   odd_re_2*sin_ang42 +
                                odd_im_3*cos_ang43     +   odd_re_3*sin_ang43 +
                                odd_im_4*cos_ang44     +   odd_re_4*sin_ang44;     

            float cos_ang_shift4 = cosf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10));
            float sin_ang_shift4 = sinf(-2.0*M_PI*(1.0/lgth)*(k+4*lgth/10));


            X_out[idx_start + 8*lgth/10]          = even_sum_real + odd_sum_real*cos_ang_shift4  - odd_sum_imag*sin_ang_shift4;
            X_out[idx_start + 8*lgth/10+1]        = even_sum_imag + odd_sum_imag*cos_ang_shift4  + odd_sum_real*sin_ang_shift4;  
            X_out[idx_start + 8*lgth/10 + lgth]   = even_sum_real - (odd_sum_real*cos_ang_shift4 - odd_sum_imag*sin_ang_shift4);
            X_out[idx_start + 8*lgth/10+1 + lgth] = even_sum_imag - (odd_sum_imag*cos_ang_shift4 + odd_sum_real*sin_ang_shift4);
            
        } 

    }
}


