#include <R.h>
#include <stdlib.h>
#include "helpers.h"

/*
 * Find pairs that can be the result of over-digestion or under-digestion.
 * Reads that start at the same position but end at a different one or
 * start at a different position but end at the same one
 */

void same_end (int *end_a, int *len_a,
               int *end_b, int *len_b,
               int *out_a, int *out_b)

{
    int i, j;
    int a_pos, b_pos;
    int x = 1;

    int **sorted_a;
    int **sorted_b;

    // create arrays 2D to be able to sort them and keep track of the original order
    sorted_a = make_list_of_pairs (*len_a, end_a);
    sorted_b = make_list_of_pairs (*len_b, end_b);

    // searching for the pairs in sorted arrays allows to do in a faster way
    qsort(sorted_a, *len_a, sizeof(sorted_a[0]), cmp_fun);
    qsort(sorted_b, *len_b, sizeof(sorted_b[0]), cmp_fun);

    // find the pairs
    j = 0;
    for (i=0; i<*len_a; i++) {
        for (; (j < *len_b && sorted_b[j][0] <= sorted_a[i][0]); j++) {
            a_pos = sorted_a[i][1];
            b_pos = sorted_b[j][1];
            if (out_b[b_pos] == 0) {
                if (sorted_a[i][0] == sorted_b[j][0]) {
                    out_a[a_pos] = x;
                    out_b[b_pos] = x;
                    x++;
                    break;
                }
            }
        }
    }

    // free the 2D arrays
    free_array_2d (sorted_a, *len_a);
    free_array_2d (sorted_b, *len_b);
}

