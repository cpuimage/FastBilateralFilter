/**
 * \file imageio.h
 * \brief Implements read_image and write_image functions
 * \author Pascal Getreuer <getreuer@cmla.ens-cachan.fr>
 *
 * Copyright (c) 2010-2013, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, or the terms of the
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef _IMAGEIO_H_
#define _IMAGEIO_H_

#include <stdio.h>
#include "basic.h"

/** \brief Limit on the maximum allowed image width or height (security). */
#define MAX_IMAGE_SIZE 10000


#ifndef DOXYGEN_SHOULD_SKIP_THIS


/* Definitions for specifying image formats */
#define IMAGEIO_U8            0x0000
#define IMAGEIO_SINGLE        0x0001
#define IMAGEIO_FLOAT         IMAGEIO_SINGLE
#define IMAGEIO_float        0x0002
#define IMAGEIO_STRIP_ALPHA   0x0010
#define IMAGEIO_BGRFLIP       0x0020
#define IMAGEIO_AFLIP         0x0040
#define IMAGEIO_GRAYSCALE     0x0080
#define IMAGEIO_GRAY          IMAGEIO_GRAYSCALE
#define IMAGEIO_PLANAR        0x0100
#define IMAGEIO_COLUMNMAJOR   0x0200
#define IMAGEIO_RGB           (IMAGEIO_STRIP_ALPHA)
#define IMAGEIO_BGR           (IMAGEIO_STRIP_ALPHA | IMAGEIO_BGRFLIP)
#define IMAGEIO_RGBA          0x0000
#define IMAGEIO_BGRA          (IMAGEIO_BGRFLIP)
#define IMAGEIO_ARGB          (IMAGEIO_AFLIP)
#define IMAGEIO_ABGR          (IMAGEIO_BGRFLIP | IMAGEIO_AFLIP)

#endif /* DOXYGEN_SHOULD_SKIP_THIS */


#ifndef _CRT_SECURE_NO_WARNINGS
/** \brief Avoid MSVC warnings on using fopen */
#define _CRT_SECURE_NO_WARNINGS
#endif

int identify_image_type(char *type, const char *filename);

void *read_image(int *width, int *height,
                 const char *filename, unsigned format);

int write_image(void *image, int width, int height,
                const char *filename, unsigned format, int quality);

#endif /* _IMAGEIO_H_ */
