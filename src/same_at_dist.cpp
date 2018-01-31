#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List equals_at_dist (S4 x, S4 y, int max_dist=5)
{
	int i;  // index for xs
    int j;  // bottom index for ys being checker
    int k;  // current ys being checked

    int a = 1;  // keep track of pairs by giving them the same number
    int d;  //dyad distance

    int min_pos;  // to keep track of what position is at the smallest distance
    int min_dist;  // keep track of the smallest distance found
    int min_found;  // to remember if at least one pair close enough has been found

    IntegerVector xs = x.slot("start");
    IntegerVector ys = y.slot("start");
    int xlen = xs.size();
    int ylen = ys.size();

    IntegerVector out_x(xlen);
    IntegerVector out_y(ylen);

    for (i=0, j=0; i < xlen; ++i) {
        min_found = 0;
        // pre-set the minimum dist to a distance that will be bigger than
        // eny new found distance
        min_dist = max_dist + 1;

        if (j > ylen) {  // all ys are done
            break;
        }

        // set the first value where the distance is in the range
        while ((d = xs[i] - ys[j]) > max_dist) {
            ++j;
        }

        // check from when distance stops being too small (j) until when it's too big
        for (k = j; ((d = xs[i] - ys[i]) >= -max_dist) && k < ylen; ++k) {
            if (out_y[k] != 0) {
                continue;
            }

            min_found = 1;  // if we get here, we've found at least one matching distance

            if (abs(d) < min_dist) {  // new distance is closer than the previous one, keep it
                min_dist = abs(d);
                min_pos = k;
            }
        }

        if (min_found) {  // we found at least one pair close enough
            out_x[i] = a;
            out_y[min_pos] = a;
            ++a;
        }
    }

    return List::create(out_x, out_y);
}
