/*
 * _ADCn_INx_AR.c
 *
 * Created on: Oct 7, 2024
 * Author: Mwangi Alex. W
 *
 * The auto-regression algorithm is a feature extraction technique that intends to determine how a value in a time series
 * data depends on its past values within a window that the designer of the algorithm specifies. The number of past values (window)
 * chosen determines the order of the algorithm. That way, the algorithm accepts time series data as its input arguments and
 * produces coefficients from the Yule-walker linear equations that reflect the dependency of a particular value within a data set with its past values still
 * within the data set. The functions developed here therefore first of all compute the autocorrelation values within the data
 * set which in this case is data from the output of a moving average DSP scheme and go ahead to compute the coefficients using the
 * autocorrelation values
 */

//INCLUDES
#include "_ADCn_INx_AR.h"

//VARIABLES
float32_t AutoCorr_1[AR_ORDER + 1]; // Buffer holding R0 to R10, e.g., autocorr_buffer[0] = R0, ..., autocorr_buffer[10] = R10
float32_t AutoCorr_2[AR_ORDER + 1];
float32_t AutoCorr_3[AR_ORDER + 1];
float32_t AutoCorr_4[AR_ORDER + 1];
float32_t AutoCorr_5[AR_ORDER + 1];
float32_t AutoCorr_6[AR_ORDER + 1];

float32_t AR_Coeffs_1[AR_ORDER]; // Buffer that holds the coefficient values. Initialized to zero before computation
float32_t AR_Coeffs_2[AR_ORDER];
float32_t AR_Coeffs_3[AR_ORDER];
float32_t AR_Coeffs_4[AR_ORDER];
float32_t AR_Coeffs_5[AR_ORDER];
float32_t AR_Coeffs_6[AR_ORDER];

ADC1_IN1_MA AR_ADC1_IN1; // Declaring an instance of Auto-regression feature
ADC1_IN2_MA AR_ADC1_IN2;
ADC2_IN3_MA AR_ADC2_IN3;
ADC2_IN4_MA AR_ADC2_IN4;
ADC3_IN1_MA AR_ADC3_IN1;
ADC3_IN2_MA AR_ADC3_IN2;


arm_matrix_instance_f32 ADC1_IN1_YW_mtx; // Creating a matrix instance for the CMSIS DSP matrix initialization function
arm_matrix_instance_f32 ADC1_IN2_YW_mtx;
arm_matrix_instance_f32 ADC2_IN3_YW_mtx;
arm_matrix_instance_f32 ADC2_IN4_YW_mtx;
arm_matrix_instance_f32 ADC3_IN1_YW_mtx;
arm_matrix_instance_f32 ADC3_IN2_YW_mtx;


arm_matrix_instance_f32 ADC1_IN1_Inv_YW_mtx; // Creating an inverse matrix instance for the CMSIS DSP matrix initialization function
arm_matrix_instance_f32 ADC1_IN2_Inv_YW_mtx;
arm_matrix_instance_f32 ADC2_IN3_Inv_YW_mtx;
arm_matrix_instance_f32 ADC2_IN4_Inv_YW_mtx;
arm_matrix_instance_f32 ADC3_IN1_Inv_YW_mtx;
arm_matrix_instance_f32 ADC3_IN2_Inv_YW_mtx;

arm_matrix_instance_f32 ADC1_IN1_AC_mtx;  // Creating an auto-correlation matrix instance for the CMSIS DSP matrix initialization function
arm_matrix_instance_f32 ADC1_IN2_AC_mtx;
arm_matrix_instance_f32 ADC2_IN3_AC_mtx;
arm_matrix_instance_f32 ADC2_IN4_AC_mtx;
arm_matrix_instance_f32 ADC3_IN1_AC_mtx;
arm_matrix_instance_f32 ADC3_IN2_AC_mtx;

arm_matrix_instance_f32 ADC1_IN1_coeffs_mtx;   // Creating a coefficients matrix instance for the CMSIS DSP matrix initialization function
arm_matrix_instance_f32 ADC1_IN2_coeffs_mtx;
arm_matrix_instance_f32 ADC2_IN3_coeffs_mtx;
arm_matrix_instance_f32 ADC2_IN4_coeffs_mtx;
arm_matrix_instance_f32 ADC3_IN1_coeffs_mtx;
arm_matrix_instance_f32 ADC3_IN2_coeffs_mtx;

arm_status StatusCoeffs_1; // Variable to monitor the status of the coefficients calculation operation
arm_status StatusCoeffs_2;
arm_status StatusCoeffs_3;
arm_status StatusCoeffs_4;
arm_status StatusCoeffs_5;
arm_status StatusCoeffs_6;

//FUNCTION DEFINITIONS
void ADC1_IN1_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n = 0; n <= AR_ORDER; n++)  // AR_ORDER is 10, so n runs from 0 to 10 (inclusive)
	{
		float32_t result_1 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag  for the data set*/
		uint32_t Blocksize_1 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC1_IN1.MA_ADC1_IN1_OutBfr[n]), AR_ADC1_IN1.MA_ADC1_IN1_OutBfr, Blocksize_1, &result_1); // A CMSIS DSP function for computing dot products

		AutoCorr_1[n] = result_1 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC1_IN1_autoreg_coeffs(void)
{
	float32_t AC_Matrix_1 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_1[n] = AutoCorr_1[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC1_IN1_AC_mtx, AR_ORDER, 1, AC_Matrix_1); // Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_1[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_1[r * 10 + c] = AutoCorr_1[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC1_IN1_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_1); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_1 [AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC1_IN1_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_1); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_1 = arm_mat_inverse_f32(&ADC1_IN1_YW_mtx, &ADC1_IN1_Inv_YW_mtx);

	if(StatusInv_1 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_1, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC1_IN1_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_1); //Initializes the coefficients Matrix

        StatusCoeffs_1 = arm_mat_mult_f32( &ADC1_IN1_Inv_YW_mtx, &ADC1_IN1_AC_mtx, &ADC1_IN1_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_1;
}



void ADC1_IN2_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n=0; n<AR_ORDER; n++)
	{
		float32_t result_2 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag */
		uint32_t Blocksize_2 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC1_IN2.MA_ADC1_IN2_OutBfr[n]), AR_ADC1_IN2.MA_ADC1_IN2_OutBfr, Blocksize_2, &result_2);

		AutoCorr_2[n] = result_2 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC1_IN2_autoreg_coeffs(void)
{
	float32_t AC_Matrix_2 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_2[n] = AutoCorr_2[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC1_IN2_AC_mtx, AR_ORDER, 1, AC_Matrix_2); //Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_2[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_2[r * 10 + c] = AutoCorr_2[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC1_IN2_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_2); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_2[AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC1_IN2_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_2); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_2 = arm_mat_inverse_f32(&ADC1_IN2_YW_mtx, &ADC1_IN2_Inv_YW_mtx);

	if(StatusInv_2 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_2, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC1_IN2_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_2); //Initializes the coefficients Matrix

        StatusCoeffs_2 = arm_mat_mult_f32( &ADC1_IN2_Inv_YW_mtx, &ADC1_IN2_AC_mtx, &ADC1_IN2_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_2;
}



void ADC2_IN3_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n=0; n<AR_ORDER; n++)
	{
		float32_t result_3 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag */
		uint32_t Blocksize_3 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC2_IN3.MA_ADC2_IN3_OutBfr[n]), AR_ADC2_IN3.MA_ADC2_IN3_OutBfr, Blocksize_3, &result_3);

		AutoCorr_3[n] = result_3 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC2_IN3_autoreg_coeffs(void)
{
	float32_t AC_Matrix_3 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_3[n] = AutoCorr_3[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC2_IN3_AC_mtx, AR_ORDER, 1, AC_Matrix_3); //Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_3[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_3[r * 10 + c] = AutoCorr_3[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC2_IN3_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_3); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_3[AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC2_IN3_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_3); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_3 = arm_mat_inverse_f32(&ADC2_IN3_YW_mtx, &ADC2_IN3_Inv_YW_mtx);

	if(StatusInv_3 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_3, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC2_IN3_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_3); //Initializes the coefficients Matrix

        StatusCoeffs_3 = arm_mat_mult_f32( &ADC2_IN3_Inv_YW_mtx, &ADC2_IN3_AC_mtx, &ADC2_IN3_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_3;
}



void ADC2_IN4_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n=0; n<AR_ORDER; n++)
	{
		float32_t result_4 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag */
		uint32_t Blocksize_4 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC2_IN4.MA_ADC2_IN4_OutBfr[n]), AR_ADC2_IN4.MA_ADC2_IN4_OutBfr, Blocksize_4, &result_4);

		AutoCorr_4[n] = result_4 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC2_IN4_autoreg_coeffs(void)
{
	float32_t AC_Matrix_4 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_4[n] = AutoCorr_4[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC2_IN4_AC_mtx, AR_ORDER, 1, AC_Matrix_4); //Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_4[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_4[r * 10 + c] = AutoCorr_4[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC2_IN4_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_4); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_4[AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC2_IN4_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_4); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_4 = arm_mat_inverse_f32(&ADC2_IN4_YW_mtx, &ADC2_IN4_Inv_YW_mtx);

	if(StatusInv_4 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_4, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC2_IN4_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_4); //Initializes the coefficients Matrix

        StatusCoeffs_4 = arm_mat_mult_f32( &ADC2_IN4_Inv_YW_mtx, &ADC2_IN4_AC_mtx, &ADC2_IN4_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_4;
}



void ADC3_IN1_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n=0; n<AR_ORDER; n++)
	{
		float32_t result_5 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag */
		uint32_t Blocksize_5 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC3_IN1.MA_ADC3_IN1_OutBfr[n]), AR_ADC3_IN1.MA_ADC3_IN1_OutBfr, Blocksize_5, &result_5);

		AutoCorr_5[n] = result_5 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC3_IN1_autoreg_coeffs(void)
{
	float32_t AC_Matrix_5 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_5[n] = AutoCorr_5[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC3_IN1_AC_mtx, AR_ORDER, 1, AC_Matrix_5); //Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_5[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_5[r * 10 + c] = AutoCorr_5[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC3_IN1_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_5); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_5[AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC3_IN1_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_5); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_5 = arm_mat_inverse_f32(&ADC3_IN1_YW_mtx, &ADC3_IN1_Inv_YW_mtx);

	if(StatusInv_5 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_1, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC3_IN1_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_5); //Initializes the coefficients Matrix

        StatusCoeffs_5 = arm_mat_mult_f32( &ADC3_IN1_Inv_YW_mtx, &ADC3_IN1_AC_mtx, &ADC3_IN1_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_5;
}



void ADC3_IN2_autocorr_calc(void)
{
	/* Iterate through each lag */
	for(uint32_t n=0; n<AR_ORDER; n++)
	{
		float32_t result_6 = 0.0f; // Initializing it to zero before calculation

		/* Number of valid points for the dot product at this lag */
		uint32_t Blocksize_6 = ADC_DMA_SIXTEENTHBUFFERSIZE - AR_ORDER;

		arm_dot_prod_f32(&(AR_ADC3_IN2.MA_ADC3_IN2_OutBfr[n]), AR_ADC3_IN2.MA_ADC3_IN2_OutBfr, Blocksize_6, &result_6);

		AutoCorr_6[n] = result_6 / ADC_DMA_SIXTEENTHBUFFERSIZE;
	}
}

float32_t* ADC3_IN2_autoreg_coeffs(void)
{
	float32_t AC_Matrix_6 [AR_ORDER];

	for(uint32_t n=0; n < AR_ORDER; n++)
	{
		AC_Matrix_6[n] = AutoCorr_6[n + 1]; // R1 to R10 corresponds to indices 1 to 10
	}

	arm_mat_init_f32(&ADC3_IN2_AC_mtx, AR_ORDER, 1, AC_Matrix_6); //Initializes the autocorrelations matrix


	float32_t Yule_Walker_Matrix_6[AR_ORDER * AR_ORDER]; // Array that hold the matrix data (100 elements) according to the Yule-Walker equations

	/* Filling the Yule-Walker matrix with the appropriate values (100 elements) from the auto-correlations buffer values */
	for (int r = 0; r < AR_ORDER; r++)
	{
	    for (int32_t c = 0; c < AR_ORDER; c++)
	    {
	        /* Access the autocorrelation buffer using the absolute difference of indices */
	    	Yule_Walker_Matrix_6[r * 10 + c] = AutoCorr_6[abs(r - c)];  // We are placing the elements in a one-dimensional array Yule-Walker matrix as if it represents a 10x10 matrix. This is done by calculating the correct index for each matrix element using (i * 10 + j)
	    }
	}

	arm_mat_init_f32(&ADC3_IN2_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_6); // Initializes the Yule-Walker matrix


	float32_t Yule_Walker_Matrix_Inv_6[AR_ORDER * AR_ORDER]; // Array that hold the inverse matrix data (100 elements)

	arm_mat_init_f32(&ADC3_IN2_Inv_YW_mtx, AR_ORDER, AR_ORDER, Yule_Walker_Matrix_Inv_6); // Initializes the inverse Yule-Walker matrix

	/* Calculate the inverse of the Yule-Walker matrix and return status of the operation */
	arm_status StatusInv_6 = arm_mat_inverse_f32(&ADC3_IN2_YW_mtx, &ADC3_IN2_Inv_YW_mtx);

	if(StatusInv_6 == ARM_MATH_SUCCESS) // Check if operation was successful
	{
		memset(AR_Coeffs_6, 0, AR_ORDER * sizeof(float32_t)); // Initializes the entire autocorrelations array to values zero

		arm_mat_init_f32(&ADC3_IN2_coeffs_mtx , AR_ORDER, 1 , AR_Coeffs_6); //Initializes the coefficients Matrix

        StatusCoeffs_6 = arm_mat_mult_f32( &ADC3_IN2_Inv_YW_mtx, &ADC3_IN2_AC_mtx, &ADC3_IN2_coeffs_mtx );
	}
	else
	{
         // Do something to indicate that the process has failed
	}

    return AR_Coeffs_6;
}
