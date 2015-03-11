#include <R.h>

/*
 * Find pairs of nucleosomes where one is contained into the other
 */

void contained (int *start_a, int *end_a, int *len_a,
                int *start_b, int *end_b, int *len_b,
                int *out_a, int *out_b)
{
    int i, j, k;

    int x = 1;
    j = 0;
    for (i=0; i<*len_a; i++) {
        // the starts are sorted, but not the ends...
        // j will be the first index of the set B to have a match on the start
        while (start_a[i] >= start_b[j]) {
            j++;
        }

        // check all the ks from the first start match until
        // the set B start is bigger than the set A start
        for (k=j; (end_a[i] > start_b[k] && k < *len_b); k++) {
            // if we already matched this one, don't touch it
            if (out_b[k] == 0) {
                if (end_a[i] > end_b[k]) { // the end also matches
                    // record the position that matches in both sets
                    out_a[i] = x;
                    out_b[k] = x;
                    x++;
                    // we already found the match, no need to look for more
                    break;
                }
            }
        }
    }
}
