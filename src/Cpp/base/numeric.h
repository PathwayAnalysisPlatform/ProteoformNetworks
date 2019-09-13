#ifndef BASE_NUMERIC_H
#define BASE_NUMERIC_H

#include <functional>
#include <utility>

namespace base {
   template<typename OI, typename S, typename T>
   constexpr void iota_n(OI wi, S n, T v) noexcept {
      for (; n != 0; --n, ++v) {
         *wi++ = v;
      }
   }

   template<typename II, typename T, typename BF, typename UF>
   constexpr T transform_reduce(II ai, II af, T v, BF&& bf, UF&& uf) noexcept {
      while (ai != af) {
         v = bf(std::move(v), uf(*ai++));
      }
      return v;
   }
}

#endif
