#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List same_start (S4 a, S4 b)
{
    int i, j;
    int x;  // keep track of the pairs by giving them the same number

    IntegerVector start_a = a.slot("start");
    IntegerVector start_b = b.slot("start");
    int len_a = start_a.size();
    int len_b = start_b.size();

    IntegerVector out_a(len_a);
    IntegerVector out_b(len_b);

    // we only need to define j once because, since the starts are sorted,
    // if one is not matched once, it won't be matched later
    x = 1;
    for (i = 0, j = 0; i < len_a; ++i) {
        for (; (start_b[j] <= start_a[i] && j < len_b); ++j) {
            if (out_b[j] == 0 && start_a[i] == start_b[j]) {
                out_a[i] = x;
                out_b[j] = x;
                ++x;
                break;  // found it, no need to keep iterating
            }
        }
    }

    return List::create(out_a, out_b);
}
