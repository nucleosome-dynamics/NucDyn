#include <R.h>
#include "helpers.h"
#include "dynprogfuns.h"

void shifts (int *xs, int *nx, int *ys, int *ny,
             int *x_left, int *y_left, int *x_right, int *y_right,
             int *max_dist, int *min_dist)
{
    int **S;

    // Dynamic programming algorithm
    S = alloc_table(*nx, *ny);
    init_table(S, *nx, *ny);

    fill_table(S, *max_dist, *min_dist,
               xs, *nx, ys, *ny);

    traceback(S, *max_dist, *min_dist,
              xs, x_right, x_left, *nx,
              ys, y_right, y_left, *ny);

    // Free some memory
    free_array_2d(S, *nx+1);

}
