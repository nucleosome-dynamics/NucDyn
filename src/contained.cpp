#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List contained (S4 a, S4 b)
{
    IntegerVector sa = a.slot("start");
    IntegerVector sb = b.slot("start");
    IntegerVector wa = a.slot("width");
    IntegerVector wb = b.slot("width");

    int na = sa.size();
    int nb = sb.size();

    IntegerVector out_a(na);
    IntegerVector out_b(nb);

    int x = 1;
    int k = 0;

    for (int i = 0; i < na; ++i) {
        // the starts are sorted, but not the ends...
        // k will be the first index of the set B to have a match on the start
        while (sa[i] > sb[k]) {
            ++k;
        }

        // check all the js from the first start match until
        // the start of the set B is bigger than the start of the set A
        for (int j=k; sa[i] + wa[i] - 1 > sb[j] && j < nb; ++j) {
            // if we already matched this one, don't touch it
            // check that the end matches too
            if (out_b[j] == 0 && sa[i] + wa[i] > sb[j] + wb[j]) {
                // record the position that matches in both sets
                out_a[i] = x;
                out_b[j] = x;
                ++x;
                // we already found it, no need to look for more
                break;
            }
        }
    }
    return List::create(out_a, out_b);
}
