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
 * @file arrayalloc.c
 * @brief Memory allocating and deallocating routines for 2D arrays of datatype float and fft_complex
 *
 * @author PRAVIN NAIR  <sreehari1390@gmail.com>
 **/


#include "headersreq.h"

float **alloc_array(int rows, int columns);

void dealloc_array_fl(float **arr, int m);

fft_complex **alloc_array_complex(int rows, int columns);

void dealloc_array_fl_complex(fft_complex **arr, int m);


/**
 * \brief Dynamically allocate 2D array of floats
 * \param rows      Number of rows
 * \param columns   Number of columns
 * \return pointer to 2D array
 *
 * This routine allocates memory in heap for a 2D
 * array of dimensions rows x columns and datatype
 * float.
 */

float **alloc_array(int rows, int columns) {
   
/* Allocate an array of pointers with size equal to number of rows */

    float **twoDary = (float **) (calloc(rows, sizeof(float *))); 

/* For each row, allocate an array with size equal to number of columns */

    for (int i = 0; i < rows; i++) {
        *(twoDary + i) = (calloc(columns, sizeof(float)));
    }
	 
    return twoDary;
}


/**
 * \brief Deallocate dynamically allocated 2D array of floats
 * \param arr       Pointer to 2D array
 * \param m         Number of rows
 *
 * This routine deallocates heap memory allocated for
 * 2D array of rows m and datatype float.
 */

void dealloc_array_fl(float **arr, int m) {
    int k;

/* Free memory corresponding to each row */

    for (k = 0; k < m; k++) {
        free(arr[k]);
    }

/* Free memory corresponding to the array of pointers to rows */

    free(arr);
}


/**
 * \brief Dynamically allocate 2D array of complex floats
 * \param rows      Number of rows
 * \param columns   Number of columns
 * \return pointer to 2D array
 *
 * This routine allocates memory in heap for a 2D
 * array of dimensions rows x columns and datatype
 * fft_complex.
 */

fft_complex **alloc_array_complex(int rows, int columns) {
 
/* Allocate an array of pointers with size equal to number of rows */

    fft_complex **twoDary = (fft_complex **) (calloc(rows, sizeof(fft_complex *)));
 
/* For each row, allocate an array with size equal to number of columns */

    for (int i = 0; i < rows; i++) {
        *(twoDary + i) = (calloc(columns, sizeof(fft_complex)));
    } 
    return twoDary;
}


/**
 * \brief Deallocate dynamically allocated 2D array of complex floats
 * \param arr       Pointer to 2D array
 * \param m         Number of rows
 *
 * This routine deallocates heap memory allocated for
 * 2D array of rows m and datatype fft_complex.
 */

void dealloc_array_fl_complex(fft_complex **arr, int m) {
    int k;

/* Free memory corresponding to each row */

    for (k = 0; k < m; k++) {
        free(arr[k]);
    }

/* Free memory corresponding to the array of pointers to rows */

    free(arr);
}

