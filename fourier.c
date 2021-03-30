#include "fourier.h"

#include <math.h>
#include <stdio.h>

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k] = 0;

        for (int j = 0; j < n; j++) {
            t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
        }
    }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    nft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
    // Se os sinals possuem tamanho = 1
    if (n == 1) {
        // não precisa multiplicar por cexp(...) porque
        // é sempre igual a 1;
        t[0] = s[0];
        return;
    }

    int k;
    double complex sp[n / 2], si[n / 2];
    double complex tp[n / 2], ti[n / 2];
    double complex c_exp;
    
    // sp[k] = s[2k] & si[k] = s[2k+1];
    for (k = 0; k < n / 2; k++) {
        sp[k] = s[2 * k];
        si[k] = s[2 * k + 1];
    }

    fft(sp, tp, n / 2, sign);
    fft(si, ti, n / 2, sign);

    for (k = 0; k < n / 2; k++) {
        c_exp = cexp(sign * 2 * PI * k * I / n);
        t[k] = tp[k] + ti[k] * c_exp;
        t[k + n / 2] = tp[k] - ti[k] * c_exp;
    }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
    fft(t, s, n, 1);

    for (int k = 0; k < n; k++) {
        s[k] /= n;
    }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    int i, j;
    double complex row[width], row_t[width];
    double complex column[height], column_t[height];

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            row[j] = matrix[i][j];
        }
        fft_forward(row, row_t, width);

        for (j = 0; j < width; j++) {
            matrix[i][j] = row_t[j];
        }
    }

    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
            column[i] = matrix[i][j];
        }
        fft_forward(column, column_t, height);

        for (i = 0; i < height; i++) {
            matrix[i][j] = column_t[i];
        }
    }
}

void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    int i, j;
    double complex row[width], row_t[width];
    double complex column_s[height], column_t[height];

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            row[j] = matrix[i][j];
        }
        fft_inverse(row, row_t, width);

        for (j = 0; j < width; j++) {
            matrix[i][j] = row_t[j];
        }
    }

    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
            column_s[i] = matrix[i][j];
        }
        fft_inverse(column_s, column_t, height);

        for (i = 0; i < height; i++) {
            matrix[i][j] = column_t[i];
        }
    }
}

void filter(double complex input[MAX_SIZE][MAX_SIZE],
            double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;

            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x] = g * input[y][x];
        }
    }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE],
               double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE],
               double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
