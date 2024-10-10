/*
 * _ADCn_INx_AR.h
 *
 *  Created on: Oct 7, 2024
 *  Author: Mwangi Alex. W
 *
 *  Header file that stores the functions declarations for the auto-regression algorithm
 */

#ifndef INC__ADCN_INX_AR_H_
#define INC__ADCN_INX_AR_H_

//INCLUDES
#include "_ADCn_INx_SSC.h"
#include <stdlib.h>

//DEFINES
#define AR_ORDER 10 // Order for the auto-regression model (how much it looks into the past)

//FUNCTION DECLARATION
void ADC1_IN1_autocorr_calc(void); // Calculates the auto-correlation between data and generates auto-correlation values that are stored in a buffer
void ADC1_IN2_autocorr_calc(void);
void ADC2_IN3_autocorr_calc(void);
void ADC2_IN4_autocorr_calc(void);
void ADC3_IN1_autocorr_calc(void);
void ADC3_IN2_autocorr_calc(void);

float32_t* ADC1_IN1_autoreg_coeffs(void); // Calculates the auto-regression coefficients using the Yule-Walker equations. The coefficients are stored in a buffer
float32_t* ADC1_IN2_autoreg_coeffs(void);
float32_t* ADC2_IN3_autoreg_coeffs(void);
float32_t* ADC2_IN4_autoreg_coeffs(void);
float32_t* ADC3_IN1_autoreg_coeffs(void);
float32_t* ADC3_IN2_autoreg_coeffs(void);

#endif /* INC__ADCN_INX_AR_H_ */
