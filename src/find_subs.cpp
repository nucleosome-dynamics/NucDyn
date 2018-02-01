#include <Rcpp.h>
#include <algorithm>
#include "helpers.h"
using namespace Rcpp;

IntegerVector get_idxs (IntegerVector x)
{
    int n = x.size();
    int non_zeros = 0;
    int i, j;

    for (i = 0; i < n; ++i) {
        if (x[i]) {
            ++non_zeros;
        }
    }
    IntegerVector out(non_zeros);
    for (i = 0, j = 0; i < n; ++i) {
        if (x[i]) {
            out[j] = x[i];
            ++j;
        }
    }
    std::sort(out.begin(), out.end());

    return out;
}

std::vector<std::vector<int>> make_pairs (IntegerVector x)
{
    int n = x.size();
    std::vector<std::vector<int>> out(n, std::vector<int>(2));
    for (int i = 0; i < n; ++i) {
        out[i][0] = x[i];
        out[i][1] = i;
    }
    return out;
}

// [[Rcpp::export]]
IntegerVector find_subs (IntegerVector subset)
{
    IntegerVector idx_array = get_idxs(subset);
    int in_len = subset.size();
    int out_len = idx_array.size();
    IntegerVector out(out_len);

    int i;
    int zeros;

    std::vector<std::vector<int>> sort_sub = make_pairs(subset);
    sort(sort_sub.begin(), sort_sub.end(), cmp_fun);

    zeros = 0;  // count the leading zeros
    for (i = 0; i < in_len; ++i) {
        if (sort_sub[i][0]) {
            break;
        } else {
            ++zeros;
        }
    }

    for (i = 0; i < out_len; ++i) {
        out[i] = sort_sub[i + zeros][1] + 1;
    }

    return out;
}
