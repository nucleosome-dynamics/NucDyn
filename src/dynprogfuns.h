int ** alloc_table (int m, int n);
void init_table (int **S, int m, int n);
void fill_table (int **S, int max_dist,
                 int *xs, int m, int *ys, int n);
void traceback (int **S, int max_dist,
                int *xs, int *x_l, int *x_r, int m,
                int *ys, int *y_l, int *y_r, int n);
double get_score (int distance, int max_dist);
