#ifndef BASE_INTRINSIC_H
#define BASE_INTRINSIC_H

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>

namespace base {
   template<typename T>
   constexpr std::size_t byte_size( ) noexcept {
      return sizeof(T);
   }

   template<typename T>
   constexpr std::size_t byte_size(const T& v) noexcept {
      return byte_size<T>( );
   }

   template<typename T>
   constexpr std::size_t bit_size( ) noexcept {
      return sizeof(T) * CHAR_BIT;
   }

   template<typename T>
   constexpr std::size_t bit_size(const T& v) noexcept {
      return bit_size<T>( );
   }


   #ifdef NTRACE
      #define TRACE(e)
   #else
      #define TRACE(e) e
   #endif

   template<bool B = false>
   constexpr void static_unreachable( ) noexcept {
      static_assert(B && false);
   }


   template<typename T>
   constexpr void construct_default(T& v) noexcept {
      new(&v) T;
   }

   template<typename T>
   constexpr void construct_value(T& v) noexcept {
      new(&v) T( );
   }

   template<typename T, typename... PARAMS>
   constexpr void construct(T& v, PARAMS&&... p) noexcept {
      new(&v) T(std::forward<PARAMS>(p)...);
   }

   template<typename T>
   constexpr void destroy(T& v) noexcept {
      v.~T( );
   }

   template<typename T>
   constexpr void raw_move(T& a, T& b) noexcept {
      destroy(a);
      construct(a, std::move(b));
   }

   template<typename T>
   constexpr void move(T& a, T& b) noexcept {
      if (&a != &b) {
         raw_move(a, b);
      }
   }

   template<typename T>
   constexpr void raw_assign(T& a, const T& b) noexcept {
      destroy(a);
      construct(a, b);
   }

   template<typename T>
   constexpr void assign(T& a, const T& b) noexcept {
      if (&a != &b) {
         raw_assign(a, b);
      }
   }

   template<typename T, typename A>
   constexpr T* dynamic_allocate(A alloc, std::size_t n) noexcept {
      return typename std::allocator_traits<A>::template rebind_alloc<T>(alloc).allocate(n);
   }

   template<typename T, typename A>
   constexpr void dynamic_deallocate(A alloc, T* p, std::size_t n) noexcept {
      typename std::allocator_traits<A>::template rebind_alloc<T>(alloc).deallocate(p, n);
   }

   template<typename T, typename A>
   constexpr T* dynamic_construct(A alloc, std::size_t n) noexcept {
      auto p = dynamic_allocate<T>(alloc, n);
      std::uninitialized_default_construct_n(p, n);
      return p;
   }

   template<typename T, typename A>
   constexpr void dynamic_destroy(A alloc, T* p, std::size_t n) noexcept {
      std::destroy_n(p, n);
      dynamic_deallocate(alloc, p, n);
   }


   template<auto N> //requires std::is_integral_v<decltype(N)>
   struct constexpr_integer {
      constexpr operator decltype(N)( ) const noexcept {
         return N;
      }
   };

   template<auto I, auto F, std::intmax_t S = (I <= F ? +1 : -1), typename T>
   constexpr void unrolled_for(T&& f) noexcept {
      if constexpr(I == F) {
         return;
      } else if constexpr(std::is_same_v<void, decltype(f(constexpr_integer<I>( )))>) {
         f(constexpr_integer<I>( ));
         unrolled_for<I + S, F, S>(f);
      } else if constexpr(std::is_same_v<bool, decltype(f(constexpr_integer<I>( )))>) {
         if (f(constexpr_integer<I>( ))) {
            unrolled_for<I + S, F, S>(f);
         }
      }
   }
}

#endif
