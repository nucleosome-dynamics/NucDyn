#include <R.h>
#include <stdlib.h>
#include <math.h>
#include "helpers.h"

/*
 * Helper functions used in the dynamic programming algorithm used by shifts.c
 */

const int GAP = -10000;
const double MAX_SCORE = 100;

double get_score (int distance, double max_dist)
{
    if (distance > max_dist) {
        return -INFINITY;
    } else {
        return MAX_SCORE - distance * MAX_SCORE / max_dist;
    }
}

int ** alloc_table (int m, int n)
{
    int **S;
    int i;

    S = malloc ((m+1) * sizeof (int *));
    for (i=0; i<m+1; i++) {
        S[i] = malloc ((n+1) * sizeof (int));
    }

    return S;
}

void init_table (int **S, int m, int n)
{
    int i;

    S[0][0] = 0;
    for (i=0; i<m+1; i++) {
        S[i][0] = i * GAP;
    }
    for (i=0; i<n+1; i++) {
        S[0][i] = i * GAP;
    }
}

void fill_table (int **S, double max_dist,
                 int *x, int m, int *x_s, int *x_e,
                 int *y, int n, int *y_s, int *y_e)
{
    int i, j;
    int a, b, c;
    double s;
    int d;
    int dyad_x, dyad_y;

    for (i=1; i<m+1; i++) {
        for (j=1; j<n+1; j++) {
            dyad_x = dyad_pos(x_s[x[i-1]], x_e[x[i-1]]);
            dyad_y = dyad_pos(y_s[y[j-1]], y_e[y[j-1]]);
            d = abs(dyad_x - dyad_y);
            s = get_score(abs(d), max_dist);

            a = S[i-1][j-1] + s;
            b = S[i-1][j] + GAP;
            c = S[i][j-1] + GAP;

            if (s != -INFINITY && a >= b && a >= c) {
                S[i][j] = a;
            } else if (b >= c) {
                S[i][j] = b;
            } else {
                S[i][j] = c;
            }
        }
    }
}

void traceback (int **S, double max_dist, int *r, int *l,
                int *x_start, int *x_end, int *x, int *x_l, int *x_r, int m,
                int *y_start, int *y_end, int *y, int *y_l, int *y_r, int n)
{
    int i, j;
    int a, b, c;
    double s;
    int d;
    int x_dyad, y_dyad;

    i = m;
    j = n;

    while (i != 0 && j != 0) {
        x_dyad = dyad_pos(x_start[x[i-1]], x_end[x[i-1]]);
        y_dyad = dyad_pos(y_start[y[j-1]], y_end[y[j-1]]);
        d = y_dyad - x_dyad;
        s = get_score(abs(d), max_dist);

        a = S[i-1][j-1] + s;
        b = S[i-1][j] + GAP;
        c = S[i][j-1] + GAP;

        if (s != -INFINITY && i > 0 && j > 0 && a == S[i][j]) {
            if (d > 0) {  // left shift
                x_l[x[i-1]] = *l;
                y_l[y[j-1]] = *l;
                (*l)++;
            } else if (d < 0) {  //right shift
                x_r[x[i-1]] = *r;
                y_r[y[j-1]] = *r;
                (*r)++;
            } else {  // I hope this never happens...
                Rprintf("an error happened during the shifts calculation\n");
                return;
            }
            i--;
            j--;
        } else if (i > 0 && b == S[i][j]) {
            i--;
        } else if (j > 0 && c == S[i][j]) {
            j--;
        }
    }
}

