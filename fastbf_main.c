/**
 * @file   fastbf_main.c
 * @brief  Main executable file.
 *
 * @author PRAVIN NAIR  <sreehari1390@gmail.com>
 *         ANMOL POPLI  <anmol.ap020@gmail.com>
 */

#include "headersreq.h"
#include "timing.h"
int main(int argc, char *argv[]) {
    int cores;
#ifdef __linux__
#ifdef _OPENMP
    cores = core_affinity();
#endif
#else
#ifdef _OPENMP
    cores = omp_get_num_procs(); /** \brief Number of Hyperthreads for MP */
#else
    cores=0;
#endif
#endif
    if (cores > 8) cores = 8;

    /* Check if there is a right call for algorithm */
    if (argc > 9) {
        printf("Too many arguments. \nSyntax is: FBF input sigmas sigmar output eps sigmaref noise_indicator sigman \n");
        return EXIT_FAILURE;
    }
    if (argc < 6) {
        printf("Not enough arguments. \nSyntax is: FBF input sigmas sigmar output eps \n");
        return EXIT_FAILURE;
    }

    program_params params; /** \brief Program parameters */

    /* Declaration parameters */
    float *input_image, *output_image;
    int columns, rows;
      double start ;
    int i, j, sigmas;
    float sigmar, eps, sigmaref;

    /* Read image, image filename input as argv[1] */
    input_image = (float *) read_image(&columns, &rows, argv[1], IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE);

    sscanf(argv[2], "%d", &sigmas); /* SigmaS input as argv[2] */
    sscanf(argv[3], "%f", &sigmar); /* SigmaR input as argv[3] */
    sscanf(argv[5], "%f", &eps); /* parameter epsilon input as argv[5] */
    if (argc >= 7)
        sscanf(argv[6], "%f", &sigmaref); /* Parameter to control stretching of the difference images as argv[6] */
    else
        sigmaref = 32;

    float **image = alloc_array(rows, columns); /** \brief Matrix to store input image */
    float **image_out = alloc_array(rows, columns); /** \brief Matrix to store output image */
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            image[i][j] = input_image[i * columns + j] * 255.0f;
        }
    }

    /* Check if noise is added */
    bool addnoise;
    if (argc > 7)
        addnoise = (bool) atof(argv[7]);
    else
        addnoise = 0;

    if (addnoise == 1) {
        float sigman;
        if (argc == 9)
            sscanf(argv[8], "%f", &sigman); /*Noise standard deviation as argv[8] */
        else {
            printf("Specify standard deviation of noise. \nSyntax is: FBF input sigmas sigmar output eps sigmaref noise_indicator sigman \n");
            return EXIT_FAILURE;
        }
        if (addgaussiannoise(image, rows, columns, sigman) != EXIT_SUCCESS) {/* adding gaussian noise */
            printf("Adding gaussian noise failed \n");
            return EXIT_FAILURE;
        }
        /* Clipping the noisy image */
        for (i = 0; i < rows; i++) {
            for (j = 0; j < columns; j++) {
                if (image[i][j] < 0.0) image[i][j] = 0.0;
                if (image[i][j] > 255.0) image[i][j] = 255.0;
            }
        }
    }
	 start = now();
    /* Shiftable Bilateral Filter applied to image and result stored in image_out */
    if (shiftableBF(rows, columns, sigmas, sigmar, image, image_out, cores, &params, eps) != EXIT_SUCCESS) {
        printf("Fast bilateral filter algorithm failed \n");
        return EXIT_FAILURE;
    }
	double time_interval = calcElapsed(start, now());
	 
    printf("Number of DFT coefficients used for approximating range kernel is %d \n", params.K);
    printf("Execution time: %f s\n", time_interval);

    /* Write image, image filename input as argv[4] */
    output_image = (float *) calloc(rows * columns, sizeof(float));
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            output_image[i * columns + j] = image_out[i][j] / 255.0f;
        }
    }
    if (write_image(output_image, columns, rows, argv[4], IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE, 100) !=
        1) {
        printf("Writing image failed \n");
        return EXIT_FAILURE;
    }
    dealloc_array_fl(image_out, rows);


    /* No noise condition */
    if (addnoise == 0) {
        /* Computing the difference image */
        float *diff_image = (float *) calloc(rows * columns, sizeof(float));
        for (i = 0; i < rows * columns; i++)
            diff_image[i] = (input_image[i] - output_image[i]);

        /* Stretching and clipping the difference image */
        float std = calculatestd(diff_image, rows * columns);
        for (i = 0; i < rows * columns; i++) {
            diff_image[i] = (128.0f / 255.0f) + (diff_image[i] * sigmaref / (std * 255));
            if (diff_image[i] < 0.0f) diff_image[i] = 0.0f;
            if (diff_image[i] > 1.0f) diff_image[i] = 1.0f;
        }

        /* Writing the difference image */
        if (write_image(diff_image, columns, rows, "difference.png",
                        IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE, 100) != 1) {
            printf("Writing image failed \n");
            return EXIT_FAILURE;
        }
        dealloc_array_fl(image, rows);
        free(diff_image);
        free(output_image);
    }
        /* Noisy conditions */
    else {

        float *noisy_image = (float *) calloc(rows * columns, sizeof(float));
        for (i = 0; i < rows; i++) {
            for (j = 0; j < columns; j++) {
                noisy_image[i * columns + j] = image[i][j] / 255.0f;
            }
        }
        dealloc_array_fl(image, rows);
        /* Generating the difference images */
        float *diff_image = (float *) calloc(rows * columns, sizeof(float));
        float *diff_noisyimage = (float *) calloc(rows * columns, sizeof(float));
        for (i = 0; i < rows * columns; i++) {
            diff_image[i] = (input_image[i] - output_image[i]);
            diff_noisyimage[i] = (noisy_image[i] - output_image[i]);
        }
        /* Stretching and clipping the difference images */
        float std = calculatestd(diff_image, rows * columns);
        float stdnoisy = calculatestd(diff_noisyimage, rows * columns);
        for (i = 0; i < rows * columns; i++) {

            diff_image[i] = (128.0f / 255.0f) + (diff_image[i] * sigmaref / (std * 255));
            if (diff_image[i] < 0.0f) diff_image[i] = 0.0;
            if (diff_image[i] > 1.0f) diff_image[i] = 1.0;

            diff_noisyimage[i] = (128.0f / 255.0f) + (diff_noisyimage[i] * sigmaref / (stdnoisy * 255));
            if (diff_noisyimage[i] < 0.0f) diff_noisyimage[i] = 0.0f;
            if (diff_noisyimage[i] > 1.0f) diff_noisyimage[i] = 1.0f;
        }
        /* Writing the noisy image */
        if (write_image(noisy_image, columns, rows, "noisy.png", IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE,
                        100) != 1) {
            printf("Writing image failed \n");
            return EXIT_FAILURE;
        }
        free(noisy_image);

        /* Writing the difference image : original - filtered */
        if (write_image(diff_image, columns, rows, "diff_noiseless.png",
                        IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE, 100) != 1) {
            printf("Writing image failed \n");
            return EXIT_FAILURE;
        }
        free(diff_image);

        /* Writing the difference image : noisy - filtered */
        if (write_image(diff_noisyimage, columns, rows, "diff_noisy.png",
                        IMAGEIO_float | IMAGEIO_PLANAR | IMAGEIO_GRAYSCALE, 100) != 1) {
            printf("Writing image failed \n");
            return EXIT_FAILURE;
        }
        free(diff_noisyimage);
        free(output_image);
    }
    free(input_image);
#ifdef __linux__

    /** \brief Arrays to store range kernels for plotting */
    float* intensity = (float*) calloc(2*params.T+1,sizeof(float));
    float* rkernel = (float*) calloc(2*params.T+1,sizeof(float));
    float* rkernel_approx = (float*) calloc(2*params.T+1,sizeof(float));
    for (i=0; i<=2*params.T; i++) {
        intensity[i] = i-params.T;
        rkernel[i]  = expf(-0.5*intensity[i]*intensity[i]/(sigmar*sigmar));
    }

    /* Plot target and approximate range kernels */
    float w0 = (2*M_PI)/(2*params.T+1);
    for (i=0; i<=2*params.T; i++) {
        rkernel_approx[i] = 0;
        for (j=0; j<params.K; j++) {
            rkernel_approx[i] += params.coeff[j]*cosf(j*w0*intensity[i]);
        }
    }
    gnuplot_ctrl *h1;
    h1 = gnuplot_init() ;
    gnuplot_cmd(h1, "set terminal png");
    gnuplot_cmd(h1,"set title \"Range Kernel Plot\"");
    gnuplot_setstyle(h1,"lines");
    gnuplot_cmd(h1, "set output \"rangekernel_plot.png\"");
    gnuplot_plot_xy(h1, intensity, rkernel, 2*params.T+1, "Target Gaussian");
    gnuplot_cmd(h1, "set output \"rangekernel_plot.png\"");
    gnuplot_plot_xy(h1, intensity, rkernel_approx, 2*params.T+1, "Fourier Approximation");
    gnuplot_close(h1);
    free(intensity);
    free(rkernel);
    free(rkernel_approx);

#endif // __linux__
    return EXIT_SUCCESS;
}
