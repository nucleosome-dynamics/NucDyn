#include<R.h>

void which_ran(int *xs, int *n, int *starts, int *ends, int *m, int *out)
{
    int i;
    int j;

    for (i=0; i<*n; i++) {
        for (j=0; j<*m; j++) {
            if (xs[i] >= starts[j] && xs[i] <= ends[j]) {
                out[i] = j + 1;
                break;
            }
        }
    }
}
