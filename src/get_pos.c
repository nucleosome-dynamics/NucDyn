#include<R.h>

int is_present(int *xs, int n, int a)
{
    int i;
    for (i=0; i<n; i++) {
        if (xs[i] == a) {
            return 0;
        }
    }
    return a;
}

int which_in_whole(int start, int end,
                   int *whole_starts, int *whole_ends, int n,
                   int *out, int nidxs)
{
    int i, a;
    for (i=0; i<n; i++) {
        if (whole_starts[i] == start && whole_ends[i] == end) {
            a = is_present(out, nidxs, i + 1);
            if (a) {
                return a;
            }
        }
    }
}

int find_idx_pos(int x, int *xs, int *starts, int *ends, int n,
                 int *whole_starts, int *whole_ends, int n_whole,
                 int *out, int nidxs)
{
    int i;
    for (i=0; i<n; i++) {
        if (x == xs[i]) {
            return which_in_whole(starts[i], ends[i],
                                  whole_starts, whole_ends, n_whole,
                                  out, nidxs);
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
                              whole_starts, whole_ends, *whole_n,
                              out, *nidxs);
    }
}
