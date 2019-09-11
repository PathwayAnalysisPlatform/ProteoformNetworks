#include "conversions.hpp"

// The result vector is sorted lexicographycally
vs convert_uss_to_vs(const uss& a_set) {
    vs result;
    result.assign(a_set.begin(), a_set.end());
    return result;
}