#include "conversions.hpp"

vs convert_uss_to_vs(const uss& a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    sort(result.begin(), result.end());
    return result;
}