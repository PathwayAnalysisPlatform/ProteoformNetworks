#ifndef BASE_INTEGER_H
#define BASE_INTEGER_H

#include "intrinsic.h"
#include <algorithm>
#include <climits>
#include <cstdint>
#include <utility>
#include <type_traits>

namespace base {
   template<typename T>
   constexpr bool is_even(const T& n) {
      return (n & 1) == 0;
   }

   template<typename T>
   constexpr bool is_odd(const T& n) {
      return (n & 1) == 1;
   }


   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto multiple_overflow(const T& n, const U& d) noexcept {
      return n % d;
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto multiple_underflow(const T& n, const U& d) noexcept {
      return (d - (n % d)) % d;
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr bool is_multiple(const T& n, const U& d) noexcept {
      return n % d == 0;
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto ceil_multiple(const T& n, const U& d) noexcept {
      return n + multiple_underflow(n, d);
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto floor_multiple(const T& n, const U& d) noexcept {
      return n - multiple_overflow(n, d);
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto next_multiple(const T& n, const U& d) noexcept {
      return ceil_multiple(n + 1, d);
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto previous_multiple(const T& n, const U& d) noexcept {
      return floor_multiple(n - 1, d);
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto floor_division(const T& n, const U& d) noexcept {
      return n / d;
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto ceil_division(const T& n, const U& d) noexcept {
      return n / d + !is_multiple(n, d);
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto full_division(const T& n, const U& d) noexcept {
      return std::pair(floor_division(n, d), multiple_underflow(n, d));
   }

   template<typename T, typename U> //requires std::is_integral_v<T> && std::is_integral_v<U>
   constexpr auto int_pow(const T& b, const U& e) noexcept {
      if (e == 0) {
         return T(1) * U(1);
      }
      auto temp = int_pow(b, e / 2);
      if (temp *= temp; e % 2 == 1) {
         temp *= b;
      }
      return temp;
   }


   constexpr std::size_t log2(unsigned char n) noexcept {
      return bit_size(n) - __builtin_clz(n) - 1;
   }

   constexpr std::size_t log2(unsigned short n) noexcept {
      return bit_size(n) - __builtin_clz(n) - 1;
   }

   constexpr std::size_t log2(unsigned int n) noexcept {
      return bit_size(n) - __builtin_clz(n) - 1;
   }

   constexpr std::size_t log2(unsigned long n) noexcept {
      return bit_size(n) - __builtin_clzl(n) - 1;
   }

   constexpr std::size_t log2(unsigned long long n) noexcept {
      return bit_size(n) - __builtin_clzll(n) - 1;
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t log2(const T& v) noexcept {
      return log2(std::make_unsigned_t<T>(v));
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr bool is_pow2(const T& n) noexcept {
      return n != 0 && (n & (n - 1)) == 0;
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T pow2(const T& n) noexcept {
      return T(1) << n;
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T ceil_pow2(const T& n) noexcept {
      return pow2(log2(n) + !is_pow2(n));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T floor_pow2(const T& n) noexcept {
      return pow2(log2(n));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T next_pow2(const T& n) noexcept {
      return ceil_pow2(n + 1);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T previous_pow2(const T& n) noexcept {
      return floor_pow2(n - 1);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bits_for_semirange(const T& n) noexcept {
      return (n != 0 ? log2(n) + !is_pow2(n) : 0);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bits_for_count(const T& n) noexcept {
      return bits_for_semirange(n + 1);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bytes_for_semirange(const T& n) noexcept {
      return ceil_division(bits_for_semirange(n), CHAR_BIT);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bytes_for_count(const T& n) noexcept {
      return ceil_division(bits_for_count(n), CHAR_BIT);
   }


   namespace impl {
      template<std::size_t B>
      struct integer_;

      template<>
      struct integer_<8> {
         using type = std::int8_t;
      };

      template<>
      struct integer_<16> {
         using type = std::int16_t;
      };

      template<>
      struct integer_<32> {
         using type = std::int32_t;
      };

      template<>
      struct integer_<64> {
         using type = std::int64_t;
      };


      template<std::size_t B>
      struct unsigned_integer_;

      template<>
      struct unsigned_integer_<8> {
         using type = std::uint8_t;
      };

      template<>
      struct unsigned_integer_<16> {
         using type = std::uint16_t;
      };

      template<>
      struct unsigned_integer_<32> {
         using type = std::uint32_t;
      };

      template<>
      struct unsigned_integer_<64> {
         using type = std::uint64_t;
      };
   }

   template<std::size_t N>
   using int_t = typename impl::integer_<N>::type;

   template<std::size_t N>
   using uint_t = typename impl::unsigned_integer_<N>::type;

   template<std::size_t N>
   using int_least_t = typename impl::integer_<std::max(ceil_pow2(std::size_t(N)), std::size_t(8))>::type;

   template<std::size_t N>
   using uint_least_t = typename impl::unsigned_integer_<std::max(ceil_pow2(std::size_t(N)), std::size_t(8))>::type;

   template<std::size_t N>
   using int_semirange_t = int_least_t<bits_for_semirange(std::size_t(N)) + 1>;

   template<std::size_t N>
   using uint_semirange_t = uint_least_t<bits_for_semirange(std::size_t(N))>;

   template<std::size_t N>
   using int_count_t = int_least_t<bits_for_count(std::size_t(N)) + 1>;

   template<std::size_t N>
   using uint_count_t = uint_least_t<bits_for_count(std::size_t(N))>;
}

#endif
