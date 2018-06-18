/*
 * Copyright (c) 2016, Pravin Nair <sreehari1390@gmail.com>
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
 * @file noisycomputations.c
 * @brief Required functions to add gaussian noise and find standard deviation of an array
 *
 * @author PRAVIN NAIR  <sreehari1390@gmail.com>
 **/

#include "headersreq.h"
#if defined(_WIN32)
#include <windows.h>
typedef DWORD pid_t;
pid_t getpid();
#define _getpid getpid
#endif

/**
 * \brief Calculating standard deviation of 1D array 
 * \param arr       1D array
 * \param length    length of array 'arr'
 * \return Standard deviation of array 'arr'
 *
 * This routine takes 1D input array 'arr' and
 * calculate standard deviation of the array
 */
float calculatestd(float *arr, int length) {
    int i;
    float mean = 0.0, std = 0.0;
    for (i = 0; i < length; i++)
        mean += arr[i];
    mean /= length;
    for (i = 0; i < length; i++)
        std += ((arr[i] - mean) * (arr[i] - mean));
    std /= length;
    std = sqrtf(std);
    return std;
}

/**
 * \brief Adding gaussian noise to input image
 * \param image     Input image
 * \param rows      Image height
 * \param columns   Image width
 * \param sigman    Standard deviation of gaussian noise to be added
 *
 * This routine adds gaussian noise of standard deviation
 * sigman to the input image
 */

int addgaussiannoise(float **image, int rows, int columns, float sigman) {
    float a, b, z;
    mt_init_genrand((unsigned long int) time(NULL) + (unsigned long int) getpid());
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            a = mt_genrand_res53();
            b = mt_genrand_res53();
            z = sigman * sqrtf(-2.0f * logf(a)) * cosf(2.0f * M_PI * b);
            image[i][j] += z;
        }
    }
    return EXIT_SUCCESS;
}



