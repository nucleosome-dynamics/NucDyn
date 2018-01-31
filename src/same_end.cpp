#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace std;

bool cmp_fun (const vector<int> a, const vector<int> b)
{
    if (a[0] == b[0]) {
        return a[1] < b[1];
    } else {
        return a[0] < b[0];
    }
}

vector<vector<int>> make_end_pairs (S4 x)
{
    IntegerVector start = x.slot("start");
    IntegerVector width = x.slot("width");
    int n = start.size();
    vector<vector<int>> out(n, vector<int>(2));
    for (int i = 0; i < n; ++i) {
        out[i][0] = start[i] + width[i] - 1;
        out[i][1] = i;
    }
    return out;
}

// [[Rcpp::export]]
List same_end (S4 a, S4 b)
{
    int i, j;
    int x = 1;

    IntegerVector a_start = a.slot("start");
    IntegerVector b_start = b.slot("start");

    int len_a = a_start.size();
    int len_b = b_start.size();

    IntegerVector out_a(len_a);
    IntegerVector out_b(len_b);

    // create arrays 2D to be able to sort them and keep track of the original order
    vector<vector<int>> as = make_end_pairs(a);
    vector<vector<int>> bs = make_end_pairs(b);

    // searching for pairs in sorted arrays allows to do it in a faster way
    sort(as.begin(), as.end(), cmp_fun);
    sort(bs.begin(), bs.end(), cmp_fun);

    // find the pairs
    for (i = 0, j = 0; i < len_a; ++i) {
        for (; j < len_b && bs[j][0] <= as[i][0]; ++j) {
            if (out_b[bs[j][1]] == 0 && as[i][0] == bs[j][0]) {
                out_a[as[i][1]] = x;
                out_b[bs[i][1]] = x;
                ++x;
                break;
            }
        }
    }

    return List::create(out_a, out_b);
}
