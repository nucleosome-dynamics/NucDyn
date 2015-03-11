#include <R.h>
#include <stdlib.h>
#include "helpers.h"
#include "dynprogfuns.h"

int * build_set_a (int *s, int *ran_start, int *ran_end, int ran_idx,
                   int *starts, int *ends, int n)
{
    int *set;
    int i;

    set = malloc(sizeof(int));
    *s = 0;

    for (i=0; i<n; i++) {
        if (starts[i] >= ran_start[ran_idx] && ends[i] <= ran_end[ran_idx]) {
            (*s)++;
            set = realloc(set, (*s) * sizeof(int));
            set[(*s)-1] = i;
        }
    }

    return set;
}

int * build_set_b (int *s, int *x_start, int *x_end, int *x_sub, int xs,
                   int *y_start, int *y_end, int ys, double max_dist)
{
    int i, j;
    int dyad_x, dyad_y;
    int *set;

    set = malloc(sizeof(int));
    *s = 0;

    for (j=0; j<ys; j++) {
        dyad_y = dyad_pos(y_start[j], y_end[j]);
        for (i=0; i<xs; i++) {
            dyad_x = dyad_pos(x_start[x_sub[i]], x_end[x_sub[i]]);
            if (abs(dyad_x - dyad_y) <= max_dist) {
                (*s)++;
                set = realloc(set, (*s) * sizeof(int));
                set[(*s)-1] = j;
                break;
            }
        }
    }

    return set;
}

void shifts (int *x_start, int *x_end, int *x_size,
             int *y_start, int *y_end, int *y_size,
             int *ran_start, int *ran_end, int *ran_size,
             int *x_left, int *y_left, int *x_right, int *y_right,
             double *max_dist)
{
    int i, j;

    int ran_idx;

    int *x_subset, *y_subset;
    int x_sub_size, y_sub_size;
    int **S;

    int r = 1;
    int l = 1;

    int shift_idx;
    int dyad_x, dyad_y;

    for (ran_idx=0; ran_idx<*ran_size; ran_idx++) {

        // Build subset x:
        // fragments from set x that are inside the current range subset
        // * also modifies x_sub_size
        x_subset = build_set_a(&x_sub_size, ran_start, ran_end, ran_idx,
                               x_start, x_end, *x_size);

        // Build subset y:
        // fragments from set y that whose dyad is at most 74 bp. far from the
        // dyad of a fragment in the subset x
        // * also modifies y_sub_size
        y_subset = build_set_b(&y_sub_size, x_start, x_end, x_subset, x_sub_size,
                               y_start, y_end, *y_size, *max_dist);

        // Dynamic programming algorithm
        S = alloc_table(x_sub_size, y_sub_size);
        init_table(S, x_sub_size, y_sub_size);

        fill_table(S, *max_dist,
                   x_subset, x_sub_size, x_start, x_end,
                   y_subset, y_sub_size, y_start, y_end);

        traceback(S, *max_dist, &r, &l,
                  x_start, x_end, x_subset, x_right, x_left, x_sub_size,
                  y_start, y_end, y_subset, y_right, y_left, y_sub_size);

        // Free some memory
        free_array_2d(S, x_sub_size+1);

        free(x_subset);
        x_subset = NULL;
        free(y_subset);
        y_subset = NULL;
    }
}
