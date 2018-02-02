#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

const int GAP = -10000;
const double MAX_SCORE = 100;

IntegerVector precise_dyad_pos (S4 x)
{
    IntegerVector start = x.slot("start");
    IntegerVector width = x.slot("width");
    int n = start.size();
    IntegerVector out(n);

    for (int i; i < n; ++i) {
        out[i] = (start[i] + (width[i] - 1)/2.0) * 10;
    }

    return out;
}

double get_score (int distance, int max_dist, int min_dist)
{
    int dist_range;
    int dist_diff;
    double scale;

    if (distance > max_dist || distance < min_dist) {
        return -INFINITY;
    } else {
        dist_range = max_dist - min_dist;
        dist_diff = distance - min_dist;
        scale = MAX_SCORE / dist_range;

        return MAX_SCORE - dist_diff * scale;
    }
}

std::vector< std::vector<int> > init_table (int nx, int ny)
{
    int i;
    std::vector< std::vector<int> > S(nx+1, std::vector<int>(ny+1));

    S[0][0] = 0;
    for (i = 0; i < nx+1; ++i) {
        S[i][0] = i * GAP;
    }
    for (i = 0; i < ny+1; ++i) {
        S[0][i] = i * GAP;
    }

    return S;
}

std::vector< std::vector<int> > fill_table (std::vector< std::vector<int> > S,
                                          IntegerVector xs, IntegerVector ys,
                                int max_dist, int min_dist)
{
    int a, b, c;
    double s;
    int d;

    int nx = xs.size();
    int ny = ys.size();

    for (int i = 1; i < nx+1; ++i) {
        for (int j = 1; j < ny+1; ++j) {
            d = abs(xs[i-1] - ys[j-1]);
            s = get_score(abs(d), max_dist, min_dist);

            a = S[i-1][j-1] + (int) round(s);
            b = S[i-1][j] + GAP;
            c = S[i][j-1] + GAP;

            if (s != -INFINITY && a >= b && a >= c) {
                S[i][j] = a;
            } else if (b >= c) {
                S[i][j] = b;
            } else {
                S[i][j] = c;
            }
        }
    }

    return S;
}

List traceback (std::vector< std::vector<int> > S,
                IntegerVector xs, IntegerVector ys,
                int max_dist, int min_dist)
{
    int nx = xs.size();
    int ny = ys.size();

    IntegerVector x_left(nx);
    IntegerVector x_right(nx);
    IntegerVector y_left(ny);
    IntegerVector y_right(ny);

    int a, b, c;
    double s;
    int d;

    int r = 1;
    int l = 1;
    int i = nx;
    int j = ny;

    while (i != 0 && j != 0) {
        d = ys[j-1] - xs[i-1];
        s = get_score(abs(d), max_dist, min_dist);

        a = S[i-1][j-1] + (int) round(s);
        b = S[i-1][j] + GAP;
        c = S[i][j-1] + GAP;

        if (s != -INFINITY && i > 0 && j > 0 && (int) round(a) == S[i][j]) {
            if (d < 0) {  // left shift
                x_left[i-1] = l;
                y_left[j-1] = l;
                ++l;
            } else if (d > 0) {  // right shift
                x_right[i-1] = r;
                y_right[j-1] = r;
                ++r;
            } else {  // I hope this never happens...
                Rprintf("an error happened during the shifts calculation\n");
                return List::create();
            }
            --i;
            --j;
        } else if (i > 0 && (int) round(b) == S[i][j]) {
            --i;
        } else if (j > 0 && (int) round(c) == S[i][j]) {
            --j;
        }
    }

    return List::create(_["left"] = List::create(x_left, y_left),
                        _["right"] = List::create(x_right, y_right));
}

// [[Rcpp::export]]
List do_shifts (S4 ra, S4 rb, int max_distance=74, int min_distance=10)
{
    int max_dist = max_distance * 10;
    int min_dist = min_distance * 10;

    IntegerVector xs = precise_dyad_pos(ra);
    IntegerVector ys = precise_dyad_pos(rb);

    std::vector< std::vector<int> > S = init_table(xs.size(), ys.size());
    S = fill_table(S, xs, ys, max_dist, min_dist);

    List res = traceback(S, xs, ys, max_dist, min_dist);

    return res;
}
