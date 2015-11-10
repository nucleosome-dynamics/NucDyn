#include <R.h>
#include "helpers.h"

void find_subs (int *subset, int *idx_array, int *in_len, int *out_len, int *out)
{
    int i;
    int zeros;  // number of zeros leading the array once sorted
    int **sort_sub;

    sort_sub = make_list_of_pairs(*in_len, subset);
    qsort(sort_sub, *in_len, sizeof(sort_sub[0]), cmp_fun);

    zeros = 0;  // count the 0s
    for (i=0; i<*in_len; i++) {
        if (sort_sub[i][0] != 0) {
            break;
        }
        zeros++;
    }

    for (i=0; i<*out_len; i++) {
        out[i] = sort_sub[i + zeros][1] + 1;
    }

    free_array_2d(sort_sub, *in_len);
}
