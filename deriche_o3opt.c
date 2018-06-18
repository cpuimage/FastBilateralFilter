/*
 * Copyright (c) 2016, Anmol Popli <anmol.ap020@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */


 /**
  * @file deriche_o3opt.c
  * @brief Fast gaussian filter using deriche IIR approximation of O(3)
  *
  * @author ANMOL POPLI <anmol.ap020@gmail.com>
  *         PRAVIN NAIR  <sreehari1390@gmail.com>
  **/

#include "headersreq.h"

void convolve_deriche2D(int rows, int columns, int sigma, fft_complex **ip_padded);

float Nc[3], Dc[3], Na[3], Da[3], scale;
int w;
//float filter[w+1];
/**
 * \brief Convolve input array with 1D Causal filter
 *        (Deriche Recursive algorithm)
 * \param in        Pointer to input array
 * \param out       Pointer to output array
 * \param datasize  Input array size
 *
 * This routine performs constant time convolution of the
 * 1D input array of complex floats with 1D Causal filter
 * of Deriche Recursive algorithm. The 1D filter is an
 * IIR filter.
 */
void convolve_dericheCausal(fft_complex *in, fft_complex *out, int datasize, float *filter) {

	int i, j;
	fft_complex zero_complex;
	zero_complex.imag = 0;
	zero_complex.real = 0;
	/* Compute first 3 output elements non-recursively */
	out[w] = zero_complex;
	for (i = 0; i < w + 1; i++) {
		out[w].real += (filter[i] * in[i].real);
		out[w].imag += (filter[i] * in[i].imag);
	}
	out[w + 1] = zero_complex;
	for (i = 0; i < w + 1; i++) {
		out[w + 1].real += (filter[i] * in[i + 1].real);
		out[w + 1].imag += (filter[i] * in[i + 1].imag);
	}
	out[w + 2] = zero_complex;
	for (i = 0; i < w + 1; i++) {
		out[w + 2].real += (filter[i] * in[i + 2].real);
		out[w + 2].imag += (filter[i] * in[i + 2].imag);
	}

	/* Recursive computation of output in forward direction using filter parameters Nc, Dc and scale */
	float invScale = 1.0f / scale;
	for (i = w + 3; i < datasize - w; i++) {
		out[i] = zero_complex;
		for (j = 0; j < 3; j++) {
			out[i].real += (((Nc[j] * in[i - (2 - j)].real)* invScale));
			out[i].imag += (((Nc[j] * in[i - (2 - j)].imag)*invScale));
			out[i].real = (out[i].real - (Dc[j] * out[i - (3 - j)].real));
			out[i].imag = (out[i].imag - (Dc[j] * out[i - (3 - j)].imag));
		}
	}

}

/**
 * \brief Convolve input array with 1D AntiCausal filter
 *        (Deriche Recursive algorithm)
 * \param in        Pointer to input array
 * \param out       Pointer to output array
 * \param datasize  Input array size
 *
 * This routine performs constant time convolution of the
 * 1D input array of complex floats with 1D AntiCausal filter
 * of Deriche Recursive algorithm. The 1D filter is an
 * IIR filter.
 */
void convolve_dericheAnticausal(fft_complex *in, fft_complex *out, int datasize, float *filter) {

	int i, j;
	fft_complex zero_complex;
	zero_complex.imag = 0;
	zero_complex.real = 0;
	/* Compute last 3 output elements non-recursively */
	out[datasize - 1 - w] = zero_complex;
	for (i = 0; i < w; i++) {
		out[datasize - 1 - w].real += (filter[i] * in[datasize - 1 - i].real);
		out[datasize - 1 - w].imag += (filter[i] * in[datasize - 1 - i].imag);
	}
	out[datasize - 2 - w] = zero_complex;
	for (i = 0; i < w; i++) {
		out[datasize - 2 - w].real += (filter[i] * in[datasize - 2 - i].real);
		out[datasize - 2 - w].imag += (filter[i] * in[datasize - 2 - i].imag);
	}
	out[datasize - 3 - w] = zero_complex;
	for (i = 0; i < w; i++) {

		out[datasize - 3 - w].real += (filter[i] * in[datasize - 3 - i].real);
		out[datasize - 3 - w].imag += (filter[i] * in[datasize - 3 - i].imag);
	}

	/* Recursive computation of output in backward direction using filter parameters Na, Da and scale */
	float invScale = 1.0f / scale;
	for (i = datasize - 4 - w; i >= w; i--) {
		out[i] = zero_complex;
		for (j = 0; j < 3; j++) {
			out[i].real += ((Na[j] * in[i + (j + 1)].real)* invScale);
			out[i].imag += ((Na[j] * in[i + (j + 1)].imag)* invScale);
			out[i].real = (out[i].real - (Da[j] * out[i + (j + 1)].real));
			out[i].imag = (out[i].imag - (Da[j] * out[i + (j + 1)].imag));
		}
	}

}

/**
 * \brief Convolve input array with 1D Gaussian filter
 *        (Deriche Recursive algorithm)
 * \param in        Pointer to input array
 * \param out       Pointer to output array
 * \param datasize  Input array size
 *
 * This routine performs constant time convolution of the
 * 1D input array of complex floats with 1D Gaussian filter
 * using Deriche Recursive algorithm. The input array is separately
 * convolved with Causal and AntiCausal filters and the results are
 * added to obtain the output array.
 */
void convolve_deriche1D(fft_complex *in, fft_complex *out, int datasize, float *filter) {
	/** \brief Array to store output of Causal filter convolution */
	fft_complex* out_causal = (fft_complex*)calloc(datasize, sizeof(fft_complex));
	convolve_dericheCausal(in, out_causal, datasize, filter);
	convolve_dericheAnticausal(in, out, datasize, filter);
 
	for (int i = 0; i < datasize; i++)
	{
		in[i].real = (out_causal[i].real + out[i].real);
		in[i].imag = (out_causal[i].imag + out[i].imag);
	}
	free(out_causal);
}

/**
 * \brief Apply 2D Gaussian filter to input image
 *        (Deriche Recursive Algorithm)
 * \param rows      Image height
 * \param columns   Image width
 * \param sigma     Gaussian kernel standard deviation
 * \param ip_padded Pointer to input image
 * \param op_padded Pointer to output image
 *
 * This routine applies 2D Gaussian filter of s.d.
 * sigma to input image ip_padded of dimensions
 * rows x columns and computes output image op_padded.
 * 1D filter is first convolved along rows and then
 * along columns. The 1D convolution is performed using
 * Deriche's fast recursive algorithm.
 */
void convolve_deriche2D(int rows, int columns, int sigma, fft_complex **ip_padded) {

	/** \brief Filter radius */
	w = 3 * sigma;
	/** \brief Array to store filter weights */
	float *filter = calloc(w + 1, sizeof(float));
	/** \brief Impulse response parameters  */
	float a0 = -0.8929f, a1 = 1.021f, b0 = 1.512f, w0 = 1.475f, c0 = 1.898f, b1 = 1.556f;

	/** \brief Transfer function coefficients */
	float n22c = c0 * expf(-2 * b0 / sigma) + (a0 * cosf(w0 / sigma) - a1 * sinf(w0 / sigma)) * expf(-(b0 + b1) / sigma);
	float n11c = -((a0 * cosf(w0 / sigma) - a1 * sinf(w0 / sigma)) * expf(-b0 / sigma) + a0 * expf(-b1 / sigma) +
		2 * c0 * expf(-b0 / sigma) * cosf(w0 / sigma));
	float n00c = a0 + c0;
	float d33c = -expf(-(2 * b0 + b1) / sigma);
	float d22c = expf(-2 * b0 / sigma) + 2 * expf(-(b0 + b1) / sigma) * cosf(w0 / sigma);
	float d11c = -(2 * expf(-b0 / sigma) * cosf(w0 / sigma) + expf(-b1 / sigma));
	float d11a = d11c, d22a = d22c, d33a = d33c;
	float n33a = -d33c * n00c;
	float n22a = n22c - d22c * n00c;
	float n11a = n11c - d11c * n00c;
	Nc[0] = n22c;
	Nc[1] = n11c;
	Nc[2] = n00c;
	Dc[0] = d33c;
	Dc[1] = d22c;
	Dc[2] = d11c;
	Na[0] = n11a;
	Na[1] = n22a;
	Na[2] = n33a;
	Da[0] = d11a;
	Da[1] = d22a;
	Da[2] = d33a;
	/** \brief Scale to normalize filter weights */
	scale = (Nc[0] + Nc[1] + Nc[2]) / (1 + Dc[0] + Dc[1] + Dc[2]) +
		(Na[0] + Na[1] + Na[2]) / (1 + Da[0] + Da[1] + Da[2]);

	/* Compute normalized filter weights */

	for (int i = 0; i < w + 1; i++) {
		float gnum = -(i - w) * (i - w);
		float gden = 2 * sigma * sigma;
		filter[i] = expf(gnum / gden) / scale;
	}

	/* Symmetric padding of input image with padding width equal to the filter radius w */
	symmetric_padding(rows, columns, ip_padded, w);

	/* Convolve each row with 1D Gaussian filter */
	fft_complex *out_t = calloc(columns + (2 * w), sizeof(fft_complex));
	for (int i = 0; i < rows + 2 * w; i++) convolve_deriche1D(ip_padded[i], out_t, columns + 2 * w, filter);
	free(out_t);
	fft_complex *intemp = calloc(rows + (2 * w), sizeof(fft_complex)), *outtemp = calloc(rows + (2 * w),
		sizeof(fft_complex));
	for (int j = w; j < columns + w; j++) {
		/* Convolve each column with 1D Gaussian filter */
		for (int i = 0; i < rows + (2 * w); i++) intemp[i] = ip_padded[i][j];
		convolve_deriche1D(intemp, outtemp, rows + 2 * w, filter);
		/* Store the convolved column in row of output matrix*/
		for (int i = 0; i < rows + (2 * w); i++) ip_padded[i][j] = intemp[i];
	}
	free(filter);
	free(intemp);
	free(outtemp);
}

