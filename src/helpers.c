#include <stdlib.h>

/*
 * Several helper functions used by differend other functions in the package
 */

/*
int dyad_pos (int start, int end)
{ // given a start and an end, return the center position of that dyad
    return (start + end)/2;
}
*/

int cmp_fun (const void *pa, const void *pb)
{ // function used to compare in the sorting process
    const int *a = *(const int **)pa;
    const int *b = *(const int **)pb;

    return a[0] - b[0];
}

int ** make_list_of_pairs (int x_size, int init_vals[])
{ // make a 2D array containing every end value and its original position
    int **array_2d;
    int i;

    array_2d = (int **) malloc (x_size * sizeof(int *));

    for (i=0; i<x_size; i++) {
        //reserve memory for the pair
        array_2d[i] = (int *) malloc (2 * sizeof(int));

        //initialize the pair
        array_2d[i][0] = init_vals[i];
        array_2d[i][1] = i;
    }

    return array_2d;
}

void free_array_2d (int **a, int s)
{
    int i;

    for (i=0; i<s; i++) {
        free(a[i]);
        a[i] = NULL;
    }
    free(a);
    a = NULL;
}
