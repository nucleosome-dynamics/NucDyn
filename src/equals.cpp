#include <Rcpp.h>
using namespace Rcpp;

/*
 * Find pairs of reads that start and end at the same position
 */

// [[Rcpp::export]]
List equals (S4 a, S4 b)
{
    IntegerVector start_a = a.slot("start");
    IntegerVector start_b = b.slot("start");
    IntegerVector width_a = a.slot("width");
    IntegerVector width_b = b.slot("width");

    int len_a = start_a.size();
    int len_b = start_b.size();
    IntegerVector out_a(len_a);
    IntegerVector out_b(len_b);

    int x = 1;
    int k = 0;
    for (int i = 0; i < len_a; ++i) {
        // the starts are sorted, but not the ends...
        // k will be the first index of the set B to have a match on the start
        while (start_a[i] > start_b[k]) {
            k++;
        }

        // check all the js from the first start match until the set B start is bigger
        for (int j = k; start_a[i] == start_b[j] && j < len_b; ++j) {
            // if we already matched this one, don't touch it
            if (out_b[j] == 0 && width_a[i] == width_b[j]) {
                // record the position that matches in both sets
                out_a[i] = x;
                out_b[j] = x;
                ++x;
                break;  // we already found the match, no need to look for more
            }
        }
    }

    return List::create(out_a, out_b);
}
