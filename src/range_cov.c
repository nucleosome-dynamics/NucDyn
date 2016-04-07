#include<R.h>

void cov_ranges(int *xs, int *start, int *end, int *n, int *m, double *out)
{
    int *pos_counts;

    int i;
    int j;

    pos_counts = (int*) calloc(*m, sizeof(int));

    for (i = 0; i < *n; i++) {
        for (j = start[i]-1; (j < end[i] & j < *m); j++) {
            out[j] += xs[i];
            pos_counts[j]++;
        }
    }

    for (i = 0; i < *m; i++) {
        if (pos_counts[i] == 0) {
            out[i] = 1;
        } else {
            out[i] /= pos_counts[i];
        }
    }
}
