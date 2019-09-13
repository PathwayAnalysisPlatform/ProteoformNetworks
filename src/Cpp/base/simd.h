#ifndef SIMD_H
#define SIMD_H

#include "bitset.h"
#include "functional.h"
#include "integer.h"
#include "numeric.h"
#include <algorithm>
#include <array>
#include <climits>
#include <type_traits>
#include <utility>
#include <x86intrin.h>

namespace base {
   template<typename T, std::size_t N> requires std::is_arithmetic_v<T>
   class vectorized_array {
   public:
      using block_type = T;

      static constexpr std::size_t size( ) noexcept {
         return N;
      }

      constexpr auto& operator[](std::size_t i) noexcept {
         return v_[i];
      }

      constexpr const auto& operator[](std::size_t i) const noexcept {
         return v_[i];
      }

      constexpr auto& get_builtin( ) noexcept {
         return v_;
      }

      constexpr const auto& get_builtin( ) const noexcept {
         return v_;
      }

      constexpr auto block_begin( ) noexcept {
         return &v_[0];
      }

      constexpr auto block_begin( ) const noexcept {
         return &v_[0];
      }

      constexpr auto block_end( ) noexcept {
         return &v_[0] + size( );
      }

      constexpr auto block_end( ) const noexcept {
         return &v_[0] + size( );
      }

      constexpr auto& to_bitset( ) noexcept {
         return reinterpret_cast<bitset<N * bit_size<block_type>( ), block_type>&>(*this);
      }

      constexpr const auto& to_bitset( ) const noexcept {
         return reinterpret_cast<const bitset<N * bit_size<block_type>( ), block_type>&>(*this);
      }

   private:
      T __attribute__((vector_size(N * sizeof(T)))) v_;
   };

   using simd8i   = vectorized_array<std::int8_t , 1>;
   using simd16i  = vectorized_array<std::int16_t, 1>;
   using simd32i  = vectorized_array<std::int32_t, 1>;
   using simd64i  = vectorized_array<decltype(  __m64{ }[0]), sizeof(  __m64) / sizeof(decltype(  __m64{ }[0]))>;
   using simd128i = vectorized_array<decltype(__m128i{ }[0]), sizeof(__m128i) / sizeof(decltype(__m128i{ }[0]))>;
   using simd256i = vectorized_array<decltype(__m256i{ }[0]), sizeof(__m256i) / sizeof(decltype(__m256i{ }[0]))>;
   using simd512i = vectorized_array<decltype(__m512i{ }[0]), sizeof(__m512i) / sizeof(decltype(__m512i{ }[0]))>;

   template<typename T>
   constexpr bool has_native_simd_v = false;

   template<>
   constexpr bool has_native_simd_v<simd8i> = true;

   template<>
   constexpr bool has_native_simd_v<simd16i> = true;

   template<>
   constexpr bool has_native_simd_v<simd32i> = true;

   #if defined(__MMX__)
      template<>
      constexpr bool has_native_simd_v<simd64i> = true;
   #endif

   #if defined(__SSE2__)
      template<>
      constexpr bool has_native_simd_v<simd128i> = true;
   #endif

   #if defined(__AVX2__)
      template<>
      constexpr bool has_native_simd_v<simd256i> = true;
   #endif

   #if defined(__AVX512F__)
      template<>
      constexpr bool has_native_simd_v<simd512i> = true;
   #endif

   namespace impl {
      template<std::size_t N>
      constexpr auto escoge_secuencia_( ) noexcept {
         casos: if constexpr(base::is_multiple(N, 512) && has_native_simd_v<simd512i>) {
            return simd512i{ };
         } else if constexpr(base::is_multiple(N, 256) && has_native_simd_v<simd256i>) {
            return simd256i{ };
         } else if constexpr(base::is_multiple(N, 128) && has_native_simd_v<simd128i>) {
            return simd128i{ };
         } else if constexpr(base::is_multiple(N, 64) && has_native_simd_v<simd64i>) {
            return simd64i{ };
         } else if constexpr(base::is_multiple(N, 32) && has_native_simd_v<simd32i>) {
            return simd32i{ };
         } else if constexpr(base::is_multiple(N, 16) && has_native_simd_v<simd16i>) {
            return simd16i{ };
         } else if constexpr(base::is_multiple(N, 8) && has_native_simd_v<simd8i>) {
            return simd8i{ };
         } else {
            static_unreachable( );
         }
      }
   }

   template<std::size_t N>
   class wide_scalar {
   public:
      using simd_type = typename std::decay_t<decltype(impl::escoge_secuencia_<N>( ))>;
      using block_type = typename simd_type::block_type;

      // BUG GCC
      constexpr wide_scalar( ) noexcept = default;

      // BUG GCC
      constexpr wide_scalar(const wide_scalar& w) noexcept {
         for (int i = 0; i < arr_.size( ); ++i) {
            arr_[i] = w.arr_[i];
         }
      }

      // BUG GCC
      constexpr wide_scalar& operator=(const wide_scalar& w) noexcept {
         for (int i = 0; i < arr_.size( ); ++i) {
            arr_[i] = w.arr_[i];
         }
         return *this;
      }

      constexpr auto& simd_sequence( ) noexcept {
         return arr_;
      }

      constexpr const auto& simd_sequence( ) const noexcept {
         return arr_;
      }

      constexpr auto block_begin( ) noexcept {
         return simd_sequence( ).front( ).block_begin( );
      }

      constexpr auto block_begin( ) const noexcept {
         return simd_sequence( ).front( ).block_begin( );
      }

      constexpr auto block_end( ) noexcept {
         return simd_sequence( ).back( ).block_end( );
      }

      constexpr auto block_end( ) const noexcept {
         return simd_sequence( ).back( ).block_end( );
      }

      constexpr auto& to_bitset( ) noexcept {
         return reinterpret_cast<bitset<N, block_type>&>(*this);
      }

      constexpr const auto& to_bitset( ) const noexcept {
         return reinterpret_cast<const bitset<N, block_type>&>(*this);
      }

   private:
      std::array<simd_type, N / base::bit_size<simd_type>( )> arr_;
   };

   template<std::size_t N>
   constexpr wide_scalar<N>& operator<<=(wide_scalar<N>& w, std::size_t i) noexcept {
      return w.to_bitset( ) <<= i, w;
   }

   template<std::size_t N>
   constexpr wide_scalar<N>& operator>>=(wide_scalar<N>& w, std::size_t i) noexcept {
      return w.to_bitset( ) >>= i, w;
   }

   template<std::size_t N>
   constexpr wide_scalar<N>& operator&=(wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      return std::transform(a.simd_sequence( ).begin( ), a.simd_sequence( ).end( ), b.simd_sequence( ).begin( ), a.simd_sequence( ).begin( ), [](const auto& s1, const auto& s2) noexcept {
         return reinterpret_cast<const typename wide_scalar<N>::simd_type&>(lift(s1.get_builtin( ) & s2.get_builtin( )));
      }), a;
   }

   template<std::size_t N>
   constexpr wide_scalar<N>& operator|=(wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      return std::transform(a.simd_sequence( ).begin( ), a.simd_sequence( ).end( ), b.simd_sequence( ).begin( ), a.simd_sequence( ).begin( ), [](auto s1, const auto& s2) noexcept {
         return reinterpret_cast<const typename wide_scalar<N>::simd_type&>(lift(s1.get_builtin( ) | s2.get_builtin( )));
      }), a;
   }

   template<std::size_t N>
   constexpr wide_scalar<N>& operator^=(wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      return std::transform(a.simd_sequence( ).begin( ), a.simd_sequence( ).end( ), b.simd_sequence( ).begin( ), a.simd_sequence( ).begin( ), [](auto s1, const auto& s2) noexcept {
         return reinterpret_cast<const typename wide_scalar<N>::simd_type&>(lift(s1.get_builtin( ) ^ s2.get_builtin( )));
      }), a;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator<<(const wide_scalar<N>& w, std::size_t i) noexcept {
      auto temp = w;
      temp <<= i;
      return temp;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator>>(const wide_scalar<N>& w, std::size_t i) noexcept {
      auto temp = w;
      temp >>= i;
      return temp;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator&(const wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      auto temp = a;
      temp &= b;
      return temp;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator|(const wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      auto temp = a;
      temp |= b;
      return temp;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator^(const wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      auto temp = a;
      temp ^= b;
      return temp;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> operator~(const wide_scalar<N>& w) noexcept {
      auto temp = w;
      return std::for_each(temp.simd_sequence( ).begin( ), temp.simd_sequence( ).end( ), [](auto& s) noexcept {
         s.get_builtin( ) = ~s.get_builtin( );
      }), temp;
   }

   template<std::size_t N>
   constexpr bool operator==(const wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      return std::equal(a.simd_sequence( ).begin( ), a.simd_sequence( ).end( ), b.simd_sequence( ).begin( ), [](const auto& s1, const auto& s2) noexcept {
         if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd64i>) {
            return _m_to_int64(s1.get_builtin( )) == _m_to_int64(s2.get_builtin( ));
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd128i>) {
            return _mm_movemask_epi8(_mm_cmpeq_epi64(s1.get_builtin( ), s2.get_builtin( ))) == 0xFFFF;
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd256i>) {
            return _mm256_movemask_epi8(_mm256_cmpeq_epi64(s1.get_builtin( ), s2.get_builtin( ))) == 0xFFFFFFFF;
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd512i>) {
            return _mm512_cmpeq_epi64_mask(s1.get_builtin( ), s2.get_builtin( )) == 0b11111111;
         } else {
            return s1.to_bitset( ) == s2.to_bitset( );
         }
      });
   }

   template<std::size_t N>
   constexpr bool operator!=(const wide_scalar<N>& a, const wide_scalar<N>& b) noexcept {
      return !(a == b);
   }


   template<typename T> requires std::is_same_v<T, wide_scalar<bit_size<T>( )>>
   struct word_type<T> {
      using type = typename T::block_type;
   };

   template<typename T> requires std::is_same_v<T, wide_scalar<bit_size<T>( )>>
   constexpr std::size_t word_granularity( ) noexcept {
      return word_granularity<typename T::block_type>( );
   }


   template<typename T> requires std::is_same_v<T, wide_scalar<bit_size<T>( )>>
   constexpr T bitmask(std::size_t i, std::size_t n) noexcept {
      T w = { };
      set_n_bits(w, i, n);
      return w;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> bitmask(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      auto w = v;
      reset_n_bits(w, 0, i);
      reset_n_bits(w, i + n, N - (i + n));
      return w;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> bitmask_from(const wide_scalar<N>& v, std::size_t i) noexcept {
      auto w = v;
      reset_n_bits(w, 0, i);
      return w;
   }

   template<std::size_t N>
   constexpr wide_scalar<N> bitmask_until(const wide_scalar<N>& v, std::size_t i) noexcept {
      auto w = v;
      reset_n_bits(w, i, N - i);
      return w;
   }


   template<std::size_t N>
   constexpr auto get_n_bits(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      return bitmask_until(v >> i, n);
   }

   template<std::size_t N>
   constexpr void and_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n, const wide_scalar<N>& s) noexcept {
      v &= ~(bitmask_until(~s, n) << i);
   }

   template<std::size_t N>
   constexpr void or_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n, const wide_scalar<N>& s) noexcept {
      v |= bitmask_until(s, n) << i;
   }

   template<std::size_t N>
   constexpr void xor_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n, const wide_scalar<N>& s) noexcept {
      v ^= bitmask_until(s, n) << i;
   }


   template<std::size_t N>
   constexpr void reset_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      v.to_bitset( ).reset_n(i, n);
   }

   template<std::size_t N>
   constexpr void set_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      v.to_bitset( ).set_n(i, n);
   }

   template<std::size_t N>
   constexpr void flip_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      v.to_bitset( ).flip_n(i, n);
   }

   template<std::size_t N>
   constexpr void write_n_bits(wide_scalar<N>& v, std::size_t i, std::size_t n, const wide_scalar<N>& s) noexcept {
      reset_n_bits(v, i, n);
      or_n_bits(v, i, n, s);
   }


   template<std::size_t N>
   constexpr bool get_bit(const wide_scalar<N>& v, std::size_t i) noexcept {
      return v.to_bitset( )[i];
   }

   template<std::size_t N>
   constexpr void and_bit(wide_scalar<N>& v, std::size_t i, bool b) noexcept {
      v.to_bitset( )[i] &= b;
   }

   template<std::size_t N>
   constexpr void or_bit(wide_scalar<N>& v, std::size_t i, bool b) noexcept {
      v.to_bitset( )[i] |= b;
   }

   template<std::size_t N>
   constexpr void xor_bit(wide_scalar<N>& v, std::size_t i, bool b) noexcept {
      v.to_bitset( )[i] ^= b;
   }


   template<std::size_t N>
   constexpr void reset_bit(wide_scalar<N>& v, std::size_t i) noexcept {
      v.to_bitset( )[i].reset( );
   }

   template<std::size_t N>
   constexpr void set_bit(wide_scalar<N>& v, std::size_t i) noexcept {
      v.to_bitset( )[i].set( );
   }

   template<std::size_t N>
   constexpr void flip_bit(wide_scalar<N>& v, std::size_t i) noexcept {
      v.to_bitset( )[i].flip( );
   }

   template<std::size_t N>
   constexpr void write_bit(wide_scalar<N>& v, std::size_t i, bool b) noexcept {
      v.to_bitset( )[i] = b;
   }


   template<std::size_t N>
   constexpr void reset_bits(wide_scalar<N>& v) noexcept {
      v = wide_scalar<N>{ };
   }

   template<std::size_t N>
   constexpr void set_bits(wide_scalar<N>& v) noexcept {
      v = ~wide_scalar<N>{ };
   }

   template<std::size_t N>
   constexpr void flip_bits(wide_scalar<N>& v) noexcept {
      v = ~v;
   }


   template<std::size_t N>
   constexpr bool none_bits(const wide_scalar<N>& v) noexcept {
      return std::all_of(v.simd_sequence( ).begin( ), v.simd_sequence( ).end( ), [](const auto& s) noexcept {
         if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd64i>) {
            return _m_to_int64(s.get_builtin( )) == 0;
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd128i>) {
            return _mm_testz_si128(s.get_builtin( ), s.get_builtin( ));
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd256i>) {
            return _mm256_testz_si256(s.get_builtin( ), s.get_builtin( ));
         } else if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd512i>) {
            return _mm512_test_epi64_mask(s.get_builtin( ), s.get_builtin( )) == 0;
         } else {
            return s.to_bitset( ).none_bits( );
         }
      });
   }

   template<std::size_t N>
   constexpr bool all_bits(const wide_scalar<N>& v) noexcept {
      return none_bits(~v);
   }

   template<std::size_t N>
   constexpr bool any_bits(const wide_scalar<N>& v) noexcept {
      return !none_bits(v);
   }

   template<std::size_t N>
   constexpr bool none_n_bits(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      return v.to_bitset( ).none_n(i, n);
   }

   template<std::size_t N>
   constexpr bool all_n_bits(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      return v.to_bitset( ).all_n(i, n);
   }

   template<std::size_t N>
   constexpr bool any_n_bits(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      return v.to_bitset( ).any_n(i, n);
   }


   template<std::size_t N>
   constexpr std::size_t popcount(const wide_scalar<N>& v) noexcept {
      return transform_reduce(v.simd_sequence( ).begin( ), v.simd_sequence( ).end( ), std::size_t(0), std::plus{ }, [](const auto& s) noexcept {
         if constexpr(std::is_same_v<typename wide_scalar<N>::simd_type, simd512i>) {
            return _mm512_reduce_add_epi64(_mm512_popcnt_epi64(s.builtin( )));
         } else {
            return s.to_bitset( ).count( );
         }
      });
   }

   template<std::size_t N>
   constexpr std::size_t popcount_n(const wide_scalar<N>& v, std::size_t i, std::size_t n) noexcept {
      return v.to_bitset( ).count_n(i, n);
   }


   template<std::size_t N>
   constexpr std::size_t find_first_set_bit(const wide_scalar<N>& v) noexcept {
      return v.to_bitset( ).find_first_set( );
   }

   template<std::size_t N>
   constexpr std::size_t find_last_set_bit(const wide_scalar<N>& v) noexcept {
      return v.to_bitset( ).find_last_set( );
   }

   template<std::size_t N>
   constexpr std::size_t find_next_set_bit(const wide_scalar<N>& v, std::size_t i) noexcept {
      return v.to_bitset( ).find_next_set(i);
   }

   template<std::size_t N>
   constexpr std::size_t find_previous_set_bit(const wide_scalar<N>& v, std::size_t i) noexcept {
      return v.to_bitset( ).find_previous_set(i);
   }

   template<std::size_t N>
   constexpr std::size_t find_ith_set_bit(const wide_scalar<N>& v, std::size_t i) noexcept {
      return v.to_bitset( ).find_ith_set(i);
   }


   template<std::size_t N, typename F>
   constexpr void visit_set_bits(const wide_scalar<N>& v, F&& f) noexcept {
      return v.to_bitset( ).visit_set(f);
   }

   template<std::size_t N, typename F>
   constexpr void visit_set_n_bits(const wide_scalar<N>& v, std::size_t i, std::size_t n, F&& f) noexcept {
      return v.to_bitset( ).visit_set_n(i, n, f);
   }
}

#endif
