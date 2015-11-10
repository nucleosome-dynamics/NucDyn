#include<R.h>


void reads_involved(int *coords, int *coord_num,
                    int *starts_a, int *ends_a, int *read_num_a,
                    int *starts_b, int *ends_b, int *read_num_b,
                    double *means)
{
    int i, j;
    int a, b;

    for (i=0; i<*coord_num; i++) {
        a = 0;
        b = 0;
        for (j=0; j<*read_num_a || j<*read_num_b; j++) {
            if (j < *read_num_a && starts_a[j] <= coords[i] && ends_a[j] >= coords[i]) {
                a++;
            }
            if (j < *read_num_b && starts_b[j] <= coords[i] && ends_b[j] >= coords[i]) {
                b++;
            }
        }
        means[i] = (a+b) / 2.0;
    }
}
