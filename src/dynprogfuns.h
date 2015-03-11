int ** alloc_table (int m, int n);
void init_table (int **S, int m, int n);
void fill_table (int **S, double max_dist,
                 int *x, int m, int *x_s, int *x_e,
                 int *y, int n, int *y_s, int *y_e);
void traceback (int **S, double max_dist, int *r, int *l,
                int *x_start, int *x_end, int *x, int *x_l, int *x_r, int m,
                int *y_start, int *y_end, int *y, int *y_l, int *y_r, int n);
double get_score (int distance, double max_dist);
