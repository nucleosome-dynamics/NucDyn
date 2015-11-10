#include <R.h>

/*
 * Find pairs that can be considered the same read if they are close enough
 * Done in the case where read size has been set to a fixed value
 */

void same_at_dist (int *max_dist,
                   int *xs, int *xlen,
                   int *ys, int *ylen,
                   int *out_x, int *out_y)
{
    int i;  // index for xs
    int j;  // bottom index for ys being checker
    int k;  // current ys being checked

    int a = 1;  // keep track of pairs by giving them the same number
    int d;  //dyad distance

    int min_pos;  // to keep track of what position is at the smallest distance
    int min_dist;  // keep track of the smallest distance found
    int min_found;  // to remember if at least one pair close enough has been found

    for (i=0; i < *xlen; i++) {
        min_found = 0;
        // pre-set the minimum dist to a distance that will be bigger than
        // any new found distance
        min_dist = *max_dist + 1;

        if (j > *ylen) {  // all ys done
            break;
        }

        // set to the first value where the distance is in the range
        while ((d = xs[i] - ys[j]) > *max_dist) {
            j++;
        }

        // check from when distance stops beeing too small (j) until when it's too big
        for (k = j; ((d = xs[i] - ys[k]) >= -(*max_dist)) && k < *ylen; k++) {
            if (out_y[k] != 0) {  // if already matched, leave it alone
                continue;
            }

            min_found = 1;  // if we get here, we've found at least one matching distance

            if (abs(d) < min_dist) {  // new distance is closer than previous one, keep it
                min_dist = abs(d);
                min_pos = k;
            }
        }

        if (min_found) {  // we found at least one pair close enough
            out_x[i] = a;
            out_y[min_pos] = a;
            a++;
        }
    }
}
