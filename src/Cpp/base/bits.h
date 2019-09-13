#ifndef BASE_BITS_H
#define BASE_BITS_H

#include "functional.h"
#include "intrinsic.h"
#include <cstddef>
#include <type_traits>

namespace base {
   /*template<typename T>
   struct word_type;

   template<typename T>
   using word_type_t = typename word_type<T>::type;


   template<typename T> requires std::is_integral_v<T>
   struct word_type<T> {
      using type = T;
   };*/
      template<typename T>
      struct word_type {
         using type = T;
      };

      template<typename T>
      using word_type_t = typename word_type<T>::type;


   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t word_granularity( ) noexcept {
      return bit_size<T>( );
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr T bitmask(std::size_t i, std::size_t n) noexcept {
      return ((T(n != bit_size<T>( )) << n) - 1) << i;
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bitmask(const T& v, std::size_t i, std::size_t n) noexcept {
      return v & bitmask<T>(i, n);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bitmask_from(const T& v, std::size_t i) noexcept {
      return v & bitmask<T>(i, bit_size(v) - i);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr T bitmask_until(const T& v, std::size_t i) noexcept {
      return v & bitmask<T>(0, i);
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr T get_n_bits(const T& v, std::size_t i, std::size_t n) noexcept {
      return bitmask_until(v >> i, n);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void and_n_bits(T& v, std::size_t i, std::size_t n, const identity_t<T>& s) noexcept {
      v &= ~(bitmask_until(~s, n) << i);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void or_n_bits(T& v, std::size_t i, std::size_t n, const identity_t<T>& s) noexcept {
      v |= bitmask_until(s, n) << i;
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void xor_n_bits(T& v, std::size_t i, std::size_t n, const identity_t<T>& s) noexcept {
      v ^= bitmask_until(s, n) << i;
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr void reset_n_bits(T& v, std::size_t i, std::size_t n) noexcept {
      and_n_bits(v, i, n, T(0));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void set_n_bits(T& v, std::size_t i, std::size_t n) noexcept {
      or_n_bits(v, i, n, ~T(0));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void flip_n_bits(T& v, std::size_t i, std::size_t n) noexcept {
      xor_n_bits(v, i, n, ~T(0));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void write_n_bits(T& v, std::size_t i, std::size_t n, const identity_t<T>& s) noexcept {
      reset_n_bits(v, i, n);
      or_n_bits(v, i, n, s);
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr bool get_bit(const T& v, std::size_t i) noexcept {
      return get_n_bits(v, i, 1);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void and_bit(T& v, std::size_t i, bool b) noexcept {
      and_n_bits(v, i, 1, b);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void or_bit(T& v, std::size_t i, bool b) noexcept {
      or_n_bits(v, i, 1, b);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void xor_bit(T& v, std::size_t i, bool b) noexcept {
      xor_n_bits(v, i, 1, b);
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr void reset_bit(T& v, std::size_t i) noexcept {
      and_bit(v, i, false);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void set_bit(T& v, std::size_t i) noexcept {
      or_bit(v, i, true);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void flip_bit(T& v, std::size_t i) noexcept {
      xor_bit(v, i, true);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void write_bit(T& v, std::size_t i, bool b) noexcept {
      reset_bit(v, i);
      or_bit(v, i, b);
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr void and_bits(T& v, const identity_t<T>& s) noexcept {
      and_n_bits(v, 0, bit_size(v), s);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void or_bits(T& v, const identity_t<T>& s) noexcept {
      or_n_bits(v, 0, bit_size(v), s);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void xor_bits(T& v, const identity_t<T>& s) noexcept {
      xor_n_bits(v, 0, bit_size(v), s);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void write_bits(T& v, const identity_t<T>& s) noexcept {
      write_n_bits(v, 0, bit_size(v), s);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void reset_bits(T& v) noexcept {
      reset_n_bits(v, 0, bit_size(v));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void set_bits(T& v) noexcept {
      set_n_bits(v, 0, bit_size(v));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void flip_bits(T& v) noexcept {
      flip_n_bits(v, 0, bit_size(v));
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr bool none_bits(const T& v) noexcept {
      return v == T(0);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr bool all_bits(const T& v) noexcept {
      return v == ~T(0);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr bool any_bits(const T& v) noexcept {
      return !none_bits(v);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr bool none_n_bits(const T& v, std::size_t i, std::size_t n) noexcept {
      return none_bits(bitmask(v, i, n));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr bool all_n_bits(const T& v, std::size_t i, std::size_t n) noexcept {
      return all_bits(v | ~bitmask<T>(i, n));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr bool any_n_bits(const T& v, std::size_t i, std::size_t n) noexcept {
      return !none_n_bits(v, i, n);
   }


   constexpr std::size_t popcount(unsigned char v) noexcept {
      return __builtin_popcount(v);
   }

   constexpr std::size_t popcount(unsigned short v) noexcept {
      return __builtin_popcount(v);
   }

   constexpr std::size_t popcount(unsigned int v) noexcept {
      return __builtin_popcount(v);
   }

   constexpr std::size_t popcount(unsigned long v) noexcept {
      return __builtin_popcountl(v);
   }

   constexpr std::size_t popcount(unsigned long long v) noexcept {
      return __builtin_popcountll(v);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t popcount(const T& v) noexcept {
      return popcount(std::make_unsigned_t<T>(v));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t popcount_n(const T& v, std::size_t i, std::size_t n) noexcept {
      return popcount(bitmask(v, i, n));
   }


   constexpr void byte_swap(std::uint8_t& v) noexcept {
      return;
   }

   constexpr void byte_swap(std::uint16_t& v) noexcept {
      v = __builtin_bswap16(v);
   }

   constexpr void byte_swap(std::uint32_t& v) noexcept {
      v = __builtin_bswap32(v);
   }

   constexpr void byte_swap(std::uint64_t& v) noexcept {
      v = __builtin_bswap64(v);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr void byte_swap(T& v) noexcept {
      byte_swap(reinterpret_cast<std::make_unsigned_t<T>&>(v));
   }


   constexpr std::size_t find_first_set_bit(unsigned char v) noexcept {
      return (v != 0 ? __builtin_ctz(v) : bit_size(v));
   }

   constexpr std::size_t find_first_set_bit(unsigned short v) noexcept {
      return (v != 0 ? __builtin_ctz(v) : bit_size(v));
   }

   constexpr std::size_t find_first_set_bit(unsigned int v) noexcept {
      return (v != 0 ? __builtin_ctz(v) : bit_size(v));
   }

   constexpr std::size_t find_first_set_bit(unsigned long v) noexcept {
      return (v != 0 ? __builtin_ctzl(v) : bit_size(v));
   }

   constexpr std::size_t find_first_set_bit(unsigned long long v) noexcept {
      return (v != 0 ? __builtin_ctzll(v) : bit_size(v));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t find_first_set_bit(const T& v) noexcept {
      return find_first_set_bit(std::make_unsigned_t<T>(v));
   }


   constexpr std::size_t find_last_set_bit(unsigned char v) noexcept {
      return (v != 0 ? bit_size(v) - __builtin_clz(v) - 1 : bit_size(v));
   }

   constexpr std::size_t find_last_set_bit(unsigned short v) noexcept {
      return (v != 0 ? bit_size(v) - __builtin_clz(v) - 1 : bit_size(v));
   }

   constexpr std::size_t find_last_set_bit(unsigned int v) noexcept {
      return (v != 0 ? bit_size(v) - __builtin_clz(v) - 1 : bit_size(v));
   }

   constexpr std::size_t find_last_set_bit(unsigned long v) noexcept {
      return (v != 0 ? bit_size(v) - __builtin_clzl(v) - 1 : bit_size(v));
   }

   constexpr std::size_t find_last_set_bit(unsigned long long v) noexcept {
      return (v != 0 ? bit_size(v) - __builtin_clzll(v) - 1 : bit_size(v));
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t find_last_set_bit(const T& v) noexcept {
      return find_last_set_bit(std::make_unsigned_t<T>(v));
   }


   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t find_next_set_bit(T v, std::size_t i) noexcept {
      reset_n_bits(v, 0, i + 1);
      return find_first_set_bit(v);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t find_previous_set_bit(T v, std::size_t i) noexcept {
      reset_n_bits(v, i, bit_size(v) - i);
      return find_last_set_bit(v);
   }

   template<typename T> //requires std::is_integral_v<T>
   constexpr std::size_t find_ith_set_bit(T v, std::size_t i) noexcept {
      std::size_t res = 0;
      for (auto bits = bit_size(v); bits > 1; bits /= 2) {
         auto cuantos = popcount(bitmask_until(v, bits / 2));
         if (cuantos <= i) {
            i -= cuantos;
            v >>= bits / 2;
            res += bits / 2;
         }
      }
      return (i == 0 && get_bit(v, 0) ? res : bit_size(v));
   }

   template<typename T, typename F> //requires std::is_integral_v<T>
   constexpr void visit_set_bits(T v, F&& f) noexcept {
      while (v != 0) {
         auto bit = find_first_set_bit(v);
         reset_bit(v, bit);
         f(bit);
      }
   }

   template<typename T, typename F> //requires std::is_integral_v<T>
   constexpr void visit_set_n_bits(T v, std::size_t i, std::size_t n, F&& f) noexcept {
      visit_set_bits(bitmask(v, i, n), f);
   }
}

#endif
