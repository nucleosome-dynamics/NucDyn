#include <vector>

inline bool cmp_fun (const std::vector<int> a, const std::vector<int> b)
{
    if (a[0] == b[0]) {
        return a[1] < b[1];
    } else {
        return a[0] < b[0];
    }
}
