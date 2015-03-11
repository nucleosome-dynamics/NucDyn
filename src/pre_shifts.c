#include<R.h>
#include<stdlib.h>
#include "helpers.h"

/*
 * Find shift to the left and to the right
 */

void find_dyad_pos(int *dyads, int *starts, int *ends, int n)
{
    int i;

    for (i=0; i<n; i++) {
        dyads[i] = dyad_pos(starts[i], ends[i]);
    }
}

void find_shifts(int **dyads_a, int **dyads_b, int n_a, int n_b,
                 int *left_a, int *left_b,
                 int *right_a, int *right_b,
                 int max_diff)
{
    int big_dist;
    int best_dist, my_dist;
    int best_pos_a, best_pos_b;

    int i, j, k;
    int x, y;

    int a_val, a_idx;
    int b_val, b_idx;

    big_dist = max_diff + 1;  // any number greater than max_diff will work just fine

    x = 1;
    y = 1;

    j = 0;
    for (i=0; i<n_a; i++) {
        a_val = dyads_a[i][0];
        a_idx = dyads_a[i][1];

        while (j < n_b && (dyads_a[i][0] - dyads_b[j][0]) > max_diff) {
            j++;
        }
        best_dist = big_dist;  // reinitialize the best position to something big
        for (k=j; k<n_b; k++) {
            b_val = dyads_b[k][0];
            b_idx = dyads_b[k][1];
            // if they already have a pair, leave them alone
            if (left_b[b_idx] != 0 || right_b[b_idx] != 0) {
                continue;
            }
            my_dist = a_val - b_val;

            if (abs(my_dist) > max_diff) {
                break;
            }

            if (abs(my_dist) < abs(best_dist)) {
                best_dist = my_dist;
                best_pos_a = a_idx;
                best_pos_b = b_idx;
            }
        }

        if (best_dist == big_dist) {  // no pair found
            continue;
        }

        if (best_dist > 0) {  // shift to the left
            left_a[best_pos_a] = x;
            left_b[best_pos_b] = x;
            x++;
        } else if (best_dist < 0) {  // shift to the right
            right_a[best_pos_a] = y;
            right_b[best_pos_b] = y;
            y++;
        } else {  // I hope this never happens...
            Rprintf("there's probably something wrong with the calculations\n");
        }
    }
}

void pre_shifts(int *start_a, int *end_a, int *len_a,
            int *start_b, int *end_b, int *len_b,
            int *left_a, int *left_b,
            int *right_a, int *right_b,
            int *max_diff)
{
    int *dyads_a, *dyads_b;  // arrays with the read center position
    int **sort_dyads_a, **sort_dyads_b;  // we'll sort by dyad position

    dyads_a = (int *) malloc (*len_a * sizeof(int));
    dyads_b = (int *) malloc (*len_b * sizeof(int));

    find_dyad_pos(dyads_a, start_a, end_a, *len_a);
    find_dyad_pos(dyads_b, start_b, end_b, *len_b);

    sort_dyads_a = make_list_of_pairs (*len_a, dyads_a);
    sort_dyads_b = make_list_of_pairs (*len_b, dyads_b);

    free(dyads_a);
    free(dyads_b);
    dyads_a = NULL;
    dyads_b = NULL;

    qsort(sort_dyads_a, *len_a, sizeof(sort_dyads_a[0]), cmp_fun);
    qsort(sort_dyads_b, *len_b, sizeof(sort_dyads_b[0]), cmp_fun);

    find_shifts(sort_dyads_a, sort_dyads_b, *len_a, *len_b,
                left_a, left_b, right_a, right_b, *max_diff);

    free_array_2d(sort_dyads_a, *len_a);
    free_array_2d(sort_dyads_b, *len_b);
}

