#ifndef BASE_ATOMIC_H
#define BASE_ATOMIC_H

#include <atomic>
#include <functional>
#include <type_traits>

namespace base {
   template<typename T>
   inline std::atomic<T>& as_atomic(T& v) noexcept requires (sizeof(T) == sizeof(std::atomic<T>) && alignof(T) == alignof(std::atomic<T>)) {
      return reinterpret_cast<std::atomic<T>&>(v);
   }


   template<typename T, typename P>
   constexpr void set_min(std::atomic<T>& destino, const T& comparar, P&& pred) noexcept {
      T actual = destino.load( );
      while (pred(comparar, actual) && !destino.compare_exchange_weak(actual, comparar)) {
         continue;
      }
   }

   template<typename T>
   constexpr void set_min(std::atomic<T>& destino, const T& comparar) noexcept {
      set_min(destino, comparar, std::less{ });
   }

   template<typename T, typename P>
   constexpr T set_max(std::atomic<T>& destino, const T& comparar, P&& pred) noexcept {
      T actual = destino.load( );
      while (pred(actual, comparar) && !destino.compare_exchange_weak(actual, comparar)) {
         continue;
      }
   }

   template<typename T>
   constexpr void set_max(std::atomic<T>& destino, const T& comparar) noexcept {
      set_max(destino, comparar, std::less{ });
   }
}

#endif
