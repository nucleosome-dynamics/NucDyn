#include <R.h>

/*
 * Find pairs that can be the result of over-digestion or under-digestion.
 * Reads that start at the same position but end at a different one or
 * start at a different position but end at the same one
 */


void same_start (int *start_a, int *len_a,
                 int *start_b, int *len_b,
                 int *out_a, int *out_b)

{
    int i, j;
    int x = 1;  // keep track of the pairs by giving them the same number

    // remove the ones that start at the same position

    // we only need to define j once because, since starts are sorted,
    // if one is not matched once, it won't get matched later
    j = 0;
    for (i=0; i<*len_a; i++) {
        for(; (start_b[j] <= start_a[i] && j < *len_b); j++) {
            if (out_b[j] == 0) {  // if already matched, leave it alone
                if (start_a[i] == start_b[j]) {
                    out_a[i] = x;
                    out_b[j] = x;
                    x++;
                    break;  // found it, no need to keep iterating
                }
            }
        }
    }
}

