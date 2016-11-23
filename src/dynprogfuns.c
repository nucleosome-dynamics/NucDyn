#include <R.h>
#include "helpers.h"

/*
 * Helper functions used in the dynamic programming algorithm used by shifts.c
 */

const int GAP = -10000;
const double MAX_SCORE = 100;

double get_score (int distance, int max_dist, int min_dist)
{
    int dist_range;
    int dist_diff;
    double scale;

    if (distance > max_dist || distance < min_dist) {
        return -INFINITY;
    } else {
        // score goes from MAX_SCORE (for distance=min_dist)
        // to 0 (for distnace=max_dist)
        dist_range = max_dist - min_dist;
        dist_diff = distance - min_dist;
        scale = MAX_SCORE / dist_range;

        return MAX_SCORE - dist_diff * scale;
    }
}

int ** alloc_table (int m, int n)
{
    int **S;
    int i;

    S = (int **) malloc ((m+1) * sizeof (int *));
    for (i=0; i<m+1; i++) {
        S[i] = (int *) malloc((n+1) * sizeof(int));
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

void fill_table (int **S, int max_dist, int min_dist,
                 int *xs, int m, int *ys, int n)
{
    int i, j;
    int a, b, c;
    double s;
    int d;

    for (i=1; i<m+1; i++) {
        for (j=1; j<n+1; j++) {
            d = abs(xs[i-1] - ys[j-1]);
            s = get_score(abs(d), max_dist, min_dist);

            a = S[i-1][j-1] + (int) round(s);
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

void traceback (int **S, int max_dist, int min_dist,
                int *xs, int *x_l, int *x_r, int m,
                int *ys, int *y_l, int *y_r, int n)
{
    int i, j;
    double a, b, c;
    double s;
    int d;

    int r, l;

    r = 1;
    l = 1;

    i = m;
    j = n;

    while (i != 0 && j != 0) {
        d = ys[j-1] - xs[i-1];
        s = get_score(abs(d), max_dist, min_dist);

        a = S[i-1][j-1] + (int) round(s);
        b = S[i-1][j] + GAP;
        c = S[i][j-1] + GAP;

        if (s != -INFINITY && i > 0 && j > 0 && (int) round(a) == S[i][j]) {
            if (d > 0) {  // left shift
                x_l[i-1] = l;
                y_l[j-1] = l;
                l++;
            } else if (d < 0) {  //right shift
                x_r[i-1] = r;
                y_r[j-1] = r;
                r++;
            } else {  // I hope this never happens...
                Rprintf("an error happened during the shifts calculation\n");
                return;
            }
            i--;
            j--;
        } else if (i > 0 && (int) round(b) == S[i][j]) {
            i--;
        } else if (j > 0 && (int) round(c) == S[i][j]) {
            j--;
        }
    }
}
