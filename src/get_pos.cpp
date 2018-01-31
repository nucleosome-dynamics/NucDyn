#include <Rcpp.h>
using namespace Rcpp;

int is_present (IntegerVector xs, int a)
{
    int n = xs.size();

    for (int i = 0; i < n; ++i) {
        if (xs[i] == a) {
            return 0;
        }
    }

    return a;
}

int which_in_whole (int start, int width, S4 range, IntegerVector out)
{
    IntegerVector whole_starts = range.slot("start");
    IntegerVector whole_widths = range.slot("width");
    int n = whole_starts.size();
    int a;

    for (int i = 0; i < n; ++i) {
        if (whole_starts[i] == start && whole_widths[i] == width) {
            a = is_present(out, i+1);
            if (a) {
                return a;
            }
        }
    }

    return 0;
}

int find_idx_pos(int x, IntegerVector xs, S4 range, S4 whole_range,
                 IntegerVector out)
{
    IntegerVector starts = range.slot("start");
    IntegerVector widths = range.slot("width");
    int n = starts.size();

    for (int i = 0; i < n; ++i) {
        if (x == xs[i]) {
            return which_in_whole(starts[i], widths[i], whole_range, out);
        }
    }

    return 0;
}

// [[Rcpp::export]]
IntegerVector find_abs_pos (IntegerVector idxs, IntegerVector xs,
                            S4 sub_range, S4 whole_range)
{
    int n = idxs.size();
    IntegerVector out(n);

    for (int i = 0; i < n; ++i) {
        out[i] = find_idx_pos(idxs[i], xs, sub_range, whole_range, out);
    }

    return out;
}
