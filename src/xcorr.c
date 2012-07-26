/*--------------------------------------------------------------------
# Filename: xcorr.c
#  Purpose: Cross correlation in time domain
#  Author: Tom Richter
#  a lot of staff (normalization) copy & pasted from
#  Obspy's xcorr.c (Hansruedi Maurer, Joachim Wassermann)
# Copyright (C) 2011 T. Richter - GNU GPL 3
#---------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>

int min(int a, int b) {
    return (a <= b) ? a : b;
}

int max(int a, int b) {
    return (a >= b) ? a : b;
}

void xcorr(float *tr1, float *tr2, double *corp, int shift, int shift_zero,
        int window, int demean, int normalize, int ndat1,
        int ndat2, int ndat1d, int ndat2d)
/*
Calculates the cross-correlation function of data arrays tr1 and tr2. We use
the following definition of cross-correlation:
corp[i] = sum_j(tr1[i+j]*tr2[j])
The data is demeaned before. The result is normalized.

input:
tr1, tr2:   data arrays with length ndat1, ndat2
shift:      maximal shift
shift_zero: before cross-correlation the first data array is right-shifted
            by that amount of samples
            (The effect is a right-shift of the correlation function.
            In the formula above tr1[i+j] is substituted by tr1[i+j-shift_zero])
demean:     0 or 1. if 1 the data is demeand
normalize:  0 or 1. if 1 the cross-correlation function is normalized
            (correlation coefficient of 1 than means perfect fit)
window:     only use values in this window for demeaning and normalization
            if 0: window = min(ndat1, ndat2)
            if >0: window = this paramter
ndat1, ndat2: length of data arrays
ndat1d, ndat2d: use this length for demeaning (see source code)
            (if 0 ndat1d = ndat2 = window)
output:
corp: cross-correlation function of length 2*shift+1

 */ {
    int a, a2, b, b2, bmin, bmax, flag = 0, ind1, ind2, ind3, ind4;
    double sum, sum1, sum2, cmax;
    float *tra1, *tra2;

    tra1 = (float *) calloc(ndat1, sizeof (float));
    if (tra1 == NULL) {
        fprintf(stderr, "\nMemory allocation error!\n");
        exit(0);
    }
    tra2 = (float *) calloc(ndat2, sizeof (float));
    if (tra2 == NULL) {
        fprintf(stderr, "\nMemory allocation error!\n");
        exit(0);
    }

    /* Set standard values */
    if (window == 0) {
        window = min(ndat1, ndat2);
    }
    if (ndat1d == 0) {
        ndat1d = window;
    }
    if (ndat2d == 0) {
        ndat2d = window;
    }

    ind1 = max(0, (ndat1 - window) / 2);
    ind2 = min(ndat1, (ndat1 + window) / 2);
    ind3 = max(0, (ndat2 - window) / 2);
    ind4 = min(ndat2, (ndat2 + window) / 2);

    /* Demean data (Zero offset) */
    if (demean > 0) {
        sum = 0;
        for (a = ind1; a < ind2; a++) {
            sum += tr1[a];
        }
        sum /= ndat1d;
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a] - (float) sum;
        }
        if (sum == 0.0)
            flag = 1;
        sum = 0;
        for (a = ind3; a < ind4; a++) {
            sum += tr2[a];
        }
        sum /= ndat2d;
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a] - (float) sum;
        }
        if (sum == 0.0)
            flag += 1;
    } else {
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a];
        }
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a];
        }
    }

    /* Normalizing the traces  (max amp = 1) */
    /*if (normalize > 0) {
        cmax = 0;
        for (a = 0; a < ndat1; a++) {
            if (fabs(tra1[a]) > cmax) {
                cmax = fabs(tra1[a]);
            }
        }
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tra1[a] / (float) cmax;
        }
        cmax = 0;

        for (a = 0; a < ndat2; a++) {
            if (fabs(tra2[a]) > cmax) {
                cmax = fabs(tra2[a]);
            }
        }
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tra2[a] / (float) cmax;
        }
    }*/

    /* xcorr ... */
    //printf("%d\n", flag);
    //printf("NEW%d\n", 0);
    //printf("%d %d %d %d\n", ndat1, ndat2, shift, window);
    //printf("%d\n", 0);
    if (flag == 0) {
        a = 0;
        a2 = -shift_zero - shift;
        if (ndat1 != ndat2) {
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2 + (ndat2 - ndat1) / 2);
                bmax = min(ndat2, -a2 + (ndat1 + ndat2) / 2);
                b2 = bmin;
                b = b2 + (ndat1 - ndat2) / 2;
                corp[a] = 0;
                //printf("%d - %d %d - %d %d\n", a2, b + a2, b + a2 + bmax - bmin, b2, bmax);
                if (bmin >= bmax) {
                    continue;
                }
                for (; b2 < bmax; b++, b2++) {
                    corp[a] += tra1[b + a2] * tra2[b2];
                }
            }
        } else { // same as above only ndat2 = ndat1, we need one variable less in the second loop
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2);
                bmax = min(ndat1, -a2 + ndat1);
                corp[a] = 0;
                if (bmin >= bmax) {
                    continue;
                }
                for (b = bmin; b < bmax; b++) {
                    corp[a] += tra1[b + a2] * tra2[b];
                }
            }
        }

        /* normalize xcorr function */
        if (normalize > 0) {
            sum1 = sum2 = 0.0;
            for (a = ind1; a < ind2; a++) {
                sum1 += (*(tra1 + a))*(*(tra1 + a));
            }
            for (a = ind3; a < ind4; a++) {
                sum2 += (*(tra2 + a))*(*(tra2 + a));
            }
            sum1 = sqrt(sum1);
            sum2 = sqrt(sum2);
            cmax = 1 / (sum1 * sum2);
            for (a = 0; a < (2 * shift + 1); a++) {
                corp[a] *= cmax;
            }
        }
        /*
        // Find maximum correlation coefficient and shift
        cmax = 0;
        shift_max = 1;
        for (a=0;a<(2*shift+1);a++)
        {
            if (fabs(corp[a]) > cmax)
            {
                cmax = fabs(corp[a]);
                shift_max = a;
            }
        }
        corp_max = corp[shift_max];
        shift_max = shift_max - shift - 1;
         */
    } else {
        for (a = 0; a < (2 * shift + 1); a++) {
            corp[a] = 0;
            //shift_max = 0;
            //corp_max = 0;
        }
    }
    free((char *) tra1);
    free((char *) tra2);
}