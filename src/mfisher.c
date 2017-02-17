#include <R.h>
#include <stdlib.h>
#include <math.h>

double fisher( long a,  long b,  long c,  long d) {
    double N, za, zb, zc, zd;
    double lp;
    long max_c, min_c;
    long mode;
    double p, z, tza, tzb, tzc, tzd;
    long tc;

    za = (double) a;
    zb = (double) b;
    zc = (double) c;
    zd = (double) d;
    N = a + b + c + d;

    mode  = (long) ((zc + zd + 1.0) * (zc + za + 1.0) / (N + 2.0));
    lp = lgamma(za + zb + 1.0) + lgamma(za + zc + 1.0) + lgamma(zc + zd + 1.0)
       + lgamma(zb + zd + 1.0) - lgamma(N + 1.0)       - lgamma(za + 1.0)
       - lgamma(zb + 1.0)      - lgamma(zc + 1.0)      - lgamma(zd + 1.0);
    max_c = a < d ? a + c : c + d;
    min_c = c < b ? 0 : c - b;
    p = 1.0, z = 1.0, tza = za, tzb = zb, tzc = zc, tzd = zd;

    if (c >= mode) {
        for (tc=c+1; tc<=max_c; tc++) {
            z *= ((tza--) * (tzd--)) / ((++tzc) * (++tzb));
            p += z;
        }
        z = 1.0;
        for (tc=c-1; tc>=mode; tc--) {
            z *= ((zc--) * (zb--)) / ((++za) * (++zd));
        }
        for (; tc>=min_c; tc--) {
            z *= ((zc--) * (zb--)) / ((++za) * (++zd));
            if (z <= 1.0) {
                p += z;
                break;
            }
        }
        for (--tc; tc>=min_c && z>0.0; tc--) {
            z *= ((zc--) * (zb--)) / ((++za) * (++zd));
            p += z;
        }
    } else {
        for (tc=c-1; tc>=min_c; tc--) {
            z *= ((tzc--) * (tzb--)) / ((++tza) * (++tzd));
            p += z;
        }
        z = 1.0;
        for (tc=c+1; tc<=mode; tc++) {
            z *= ((za--) * (zd--)) / ((++zc) * (++zb));
        }
        for (; tc<=max_c; tc++) {
            z *= ((za--) * (zd--)) / ((++zc) * (++zb));
            if(z <= 1.0) {
                p += z;
                break;
            }
        }
        for (++tc; tc<=max_c && z>0.0; tc++) {
            z *= ((za--) * (zd--)) / ((++zc) * (++zb));
            p += z;
        }
    }
    lp += log(p);
    return lp < 0.0 ? lp : 0.0;
}

void get_pvals (int *x, int *y, int *n, double *out) {
    int X, Y;
    int i;
    int a, b, c, d;

    X = 0;
    Y = 0;
    for (i=0; i<*n; i++) {
        X += x[i];
        Y += y[i];
    }

    for (i=0; i<*n; i++) {
        a = X-x[i];
        b = Y-y[i];
        c = x[i];
        d = y[i];
        out[i] = exp(fisher(a, b, c, d));
    }
}
