#include<R.h>

int which_in_whole(int start, int end,
                   int *whole_starts, int *whole_ends, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (whole_starts[i] == start && whole_ends[i] == end) {
            return i + 1;
        }
    }
}

int find_idx_pos(int x, int *xs, int *starts, int *ends, int n,
                 int *whole_starts, int *whole_ends, int n_whole)
{
    int i;

    for (i=0; i<n; i++) {
        if (x == xs[i]) {
            return which_in_whole(starts[i], ends[i],
                                  whole_starts, whole_ends, n_whole);
        }
    }
}

void find_abs_pos(int *idxs, int *nidxs,
                  int *xs, int *starts, int *ends, int *n,
                  int *whole_starts, int *whole_ends, int *whole_n,
                  int *out)
{
    int i;

    for (i=0; i<*nidxs; i++) {
        out[i] = find_idx_pos(idxs[i], xs, starts, ends, *n,
                              whole_starts, whole_ends, *whole_n);
    }
}
