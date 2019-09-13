#ifndef CODING_H
#define CODING_H

#include "bits.h"
#include "intrinsic.h"

namespace base {
   constexpr std::size_t gray_code(std::size_t n) {
      return n ^ (n >> 1);
   }

   constexpr std::size_t gray_code_index(std::size_t n) {
      for (auto m = n >> 1; m != 0; m >>= 1) {
         n ^= m;
      }
      return n;
   }

   constexpr std::size_t gray_code_next_flip(std::size_t n) {
      return find_first_set_bit(gray_code(n) ^ gray_code(n + 1));
   }
}

#endif
