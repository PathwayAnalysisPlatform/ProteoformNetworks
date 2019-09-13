#ifndef BASE_BITSET_H
#define BASE_BITSET_H

#include "algorithm.h"
#include "bits.h"
#include "functional.h"
#include "integer.h"
#include "iterator.h"
#include "numeric.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

namespace base {
   namespace impl {
      template<typename T>
      class referencia_ {
      public:
         constexpr referencia_(T& m, std::size_t i) noexcept
         : memoria_(m), indice_(i) {
         }

         constexpr void operator=(bool b) noexcept {
            write_bit(memoria_, indice_, b);
         }

         constexpr void operator&=(bool b) noexcept {
            and_bit(memoria_, indice_, b);
         }

         constexpr void operator|=(bool b) noexcept {
            or_bit(memoria_, indice_, b);
         }

         constexpr void operator^=(bool b) noexcept {
            xor_bit(memoria_, indice_, b);
         }

         constexpr operator bool( ) const noexcept {
            return get( );
         }

         constexpr bool get( ) const noexcept {
            return get_bit(memoria_, indice_);
         }

         constexpr void reset( ) noexcept {
            reset_bit(memoria_, indice_);
         }

         constexpr void set( ) noexcept {
            set_bit(memoria_, indice_);
         }

         constexpr void flip( ) noexcept {
            flip_bit(memoria_, indice_);
         }

         constexpr T& block( ) noexcept {
            return memoria_;
         }

         constexpr const T& block( ) const noexcept {
            return memoria_;
         }

         constexpr const std::size_t& offset( ) const noexcept {
            return indice_;
         }

      private:
         T& memoria_;
         std::size_t indice_;
      };

      template<typename T>
      class referencia_const_ {
      public:
         constexpr referencia_const_(const T& m, std::size_t i) noexcept
         : memoria_(m), indice_(i) {
         }

         constexpr operator bool( ) const noexcept {
            return get( );
         }

         constexpr bool get( ) const noexcept {
            return get_bit(memoria_, indice_);
         }

         constexpr const T& block( ) const noexcept {
            return memoria_;
         }

         constexpr const std::size_t& offset( ) const noexcept {
            return indice_;
         }

      private:
         const T& memoria_;
         std::size_t indice_;
      };

      template<typename B>
      class bitset_base_ : public B {
      public:
         using block_type = typename B::block_type;
         using B::size;
         using B::blocks;
         using B::bits_per_block;
         using B::has_trailing_bits;
         using B::trailing_bits;
         using B::block_data;

         using reference = referencia_<block_type>;
         using const_reference = referencia_const_<block_type>;

         template<typename... PARAMS>
         constexpr bitset_base_(PARAMS&&... p) noexcept
         : B(std::forward<PARAMS>(p)...) {
         }

         constexpr bool has_resource( ) const noexcept {
            return block_begin( ) != nullptr;
         }

         constexpr void release( ) noexcept {
            return B::release( );
         }

         constexpr block_type* block_begin( ) noexcept {
            return block_data( );
         }

         constexpr block_type* block_end( ) noexcept {
            return block_data( ) + blocks( );
         }

         constexpr std::reverse_iterator<block_type*> block_rbegin( ) noexcept {
            return std::reverse_iterator(block_end( ));
         }

         constexpr std::reverse_iterator<block_type*> block_rend( ) noexcept {
            return std::reverse_iterator(block_begin( ));
         }

         constexpr const block_type* block_begin( ) const noexcept {
            return block_data( );
         }

         constexpr const block_type* block_end( ) const noexcept {
            return block_data( ) + blocks( );
         }

         constexpr std::reverse_iterator<const block_type*> block_rbegin( ) const noexcept {
            return std::reverse_iterator(block_end( ));
         }

         constexpr std::reverse_iterator<const block_type*> block_rend( ) const noexcept {
            return std::reverse_iterator(block_begin( ));
         }

         constexpr block_type get_n(std::size_t i, std::size_t n) const noexcept {
            block_type res = { }; std::size_t x = 0;
            aplica_n_(block_begin( ), i, n, [&](const auto& s) {
               or_n_bits(res, x, bits_per_block( ), s);
               x += bits_per_block( );
            }, [&](const auto& s, std::size_t i, std::size_t n) {
               or_n_bits(res, x, n, get_n_bits(s, i, n));
               x += n;
            });
            return res;
         }

         constexpr void reset( ) noexcept {
            reset_n(0, size( ));
         }

         constexpr void reset_n(std::size_t i, std::size_t n) noexcept {
            aplica_n_(block_begin( ), i, n, FUNCTOR(reset_bits), FUNCTOR(reset_n_bits));
         }

         constexpr void set( ) noexcept {
            set_n(0, size( ));
         }

         constexpr void set_n(std::size_t i, std::size_t n) noexcept {
            aplica_n_(block_begin( ), i, n, FUNCTOR(set_bits), FUNCTOR(set_n_bits));
         }

         constexpr void flip( ) noexcept {
            flip_n(0, size( ));
         }

         constexpr void flip_n(std::size_t i, std::size_t n) noexcept {
            aplica_n_(block_begin( ), i, n, FUNCTOR(flip_bits), FUNCTOR(flip_n_bits));
         }

         constexpr bool none( ) const noexcept {
            return std::all_of(block_begin( ), block_end( ), FUNCTOR(none_bits));
         }

         constexpr bool none_n(std::size_t i, std::size_t n) const noexcept {
            return examina_n_(block_begin( ), i, n, FUNCTOR(none_bits), FUNCTOR(none_n_bits));
         }

         constexpr bool all( ) const noexcept {
            return std::all_of(block_begin( ), block_end( ) - has_trailing_bits( ), FUNCTOR(all_bits)) && (!has_trailing_bits( ) || all_n_bits(*(block_end( ) - 1), 0, bits_per_block( ) - trailing_bits( )));
         }

         constexpr bool all_n(std::size_t i, std::size_t n) const noexcept {
            return examina_n_(block_begin( ), i, n, FUNCTOR(all_bits), FUNCTOR(all_n_bits));
         }

         constexpr bool any( ) const noexcept {
            return !none( );
         }

         constexpr bool any_n(std::size_t i, std::size_t n) const noexcept {
            return !none_n(i, n);
         }

         constexpr std::size_t count( ) const noexcept {
            return transform_reduce(block_begin( ), block_end( ), std::size_t(0), std::plus{ }, FUNCTOR(popcount));
         }

         constexpr std::size_t count_n(std::size_t i, std::size_t n) const noexcept {
            std::size_t res = 0;
            aplica_n_(block_begin( ), i, n, [&](const auto& s) {
               res += popcount(s);
            }, [&](const auto& s, std::size_t i, std::size_t n) {
               res += popcount_n(s, i, n);
            });
            return res;
         }

         constexpr reference operator[](std::size_t i) noexcept {
            return reference(block_begin( )[i / bits_per_block( )], i % bits_per_block( ));
         }

         constexpr const_reference operator[](std::size_t i) const noexcept {
            return const_reference(block_begin( )[i / bits_per_block( )], i % bits_per_block( ));
         }

         constexpr std::size_t find_first_set( ) const noexcept {
            return primero_desde_(block_begin( ));
         }

         constexpr std::size_t find_last_set( ) const noexcept {
            return ultimo_hasta_(block_end( ));
         }

         constexpr std::size_t find_next_set(std::size_t i) const noexcept {
            auto ai = block_begin( ) + i / bits_per_block( );
            auto rel = find_next_set_bit(*ai, i % bits_per_block( ));
            return (rel != bits_per_block( ) ? (ai - block_begin( )) * bits_per_block( ) + rel : primero_desde_(ai + 1));
         }

         constexpr std::size_t find_previous_set(std::size_t i) const noexcept {
            auto ai = block_begin( ) + i / bits_per_block( );
            auto rel = find_previous_set_bit(*ai, i % bits_per_block( ));
            return (rel != bits_per_block( ) ? (ai - block_begin( )) * bits_per_block( ) + rel : ultimo_hasta_(ai));
         }

         constexpr std::size_t find_ith_set(std::size_t i) const noexcept {
            std::size_t atras = 0;
            for (const auto& actual : range(block_begin( ), block_end( ))) {
               auto cuantos = popcount(actual);
               if (i < cuantos) {
                  return find_ith_set_bit(actual, i) + atras;
               }
               i -= cuantos, atras += bits_per_block( );
            }
            return size( );
         }

         template<typename F>
         constexpr void visit_set(F&& f) const noexcept {
            for (std::size_t i = 0; i < blocks( ); ++i) {
               visit_set_bits(block_begin( )[i], [&](auto j) {
                  f(i * bits_per_block( ) + j);
               });
            }
         }

         template<typename F>
         constexpr void visit_set_n(std::size_t i, std::size_t n, F&& f) const noexcept {
            aplica_n_(block_begin( ), i, n, [&](const auto& s) {
               visit_set_bits(s, [&](auto j) {
                  f((&s - block_begin( )) * bits_per_block( ) + j);
               });
            }, [&](const auto& s, std::size_t i, std::size_t n) {
               visit_set_n_bits(s, i, n, [&](auto j) {
                  f((&s - block_begin( )) * bits_per_block( ) + j);
               });
            });
         }

      private:
         constexpr std::size_t primero_desde_(const block_type* ini) const noexcept {
            auto iter = std::find_if(ini, block_end( ), FUNCTOR(any_bits));
            return (iter != block_end( ) ? (iter - block_begin( )) * bits_per_block( ) + find_first_set_bit(*iter) : size( ));
         }

         constexpr std::size_t ultimo_hasta_(const block_type* fin) const noexcept {
            auto iter = find_last_if(block_begin( ), fin, FUNCTOR(any_bits));
            return (iter != fin ? (iter - block_begin( )) * bits_per_block( ) + find_last_set_bit(*iter) : size( ));
         }

         template<typename BT, typename FT, typename FP> //requires std::is_same_v<BT, block_type*> || std::is_same_v<BT, const block_type*>
         static constexpr void aplica_n_(BT ini, std::size_t i, std::size_t n, FT&& ft, FP&& fp) noexcept {
            examina_n_(ini, i, n, [&](auto&&... p) {
               return ft(p...), true;
            }, [&](auto&&... p) {
               return fp(p...), true;
            });
         }

         template<typename BT, typename FT, typename FP> //requires std::is_same_v<BT, block_type*> || std::is_same_v<BT, const block_type*>
         static constexpr bool examina_n_(BT ini, std::size_t i, std::size_t n, FT&& ft, FP&& fp) noexcept {
            auto ai = ini + i / bits_per_block( );
            if (!is_multiple(i, bits_per_block( ))) {
               auto m = std::min(n, multiple_underflow(i, bits_per_block( )));
               if (!fp(*ai++, i % bits_per_block( ), m)) {
                  return false;
               }; n -= m;
            }

            auto completos = n / bits_per_block( );
            for (auto& actual : range(ai, ai + completos)) {
               if (!ft(actual)) {
                  return false;
               }
            }; n -= completos * bits_per_block( );

            if (n != 0 && !fp(ai[completos], 0, n)) {
               return false;
            }

            return true;
         }
      };

      template<typename B>
      constexpr bitset_base_<B>& operator<<=(bitset_base_<B>& b, std::size_t n) noexcept {
         if (n == b.size( )) {
            b.reset( );
         } else if (n != 0) {
            const auto wshift = n / b.bits_per_block( );
            const auto offset = n % b.bits_per_block( );
            if (offset == 0) {
               for (std::size_t i = b.blocks( ) - 1; i >= wshift; --i) {
                  b.block_begin( )[i] = b.block_begin( )[i - wshift];
               }
            } else {
               const std::size_t sub_offset = b.bits_per_block( ) - offset;
               for (std::size_t i = b.blocks( ) - 1; i > wshift; --i) {
                  b.block_begin( )[i] = (b.block_begin( )[i - wshift] << offset) | (b.block_begin( )[i - wshift - 1] >> sub_offset);
               }
               b.block_begin( )[wshift] = b.block_begin( )[0] << offset;
            }
            std::fill(b.block_begin( ), b.block_begin( ) + wshift, typename bitset_base_<B>::block_type(0));
            base::reset_n_bits(b.block_begin( )[b.blocks( ) - 1], b.bits_per_block( ) - b.trailing_bits( ), b.trailing_bits( ));
         }
         return b;
      }

      template<typename B>
      constexpr bitset_base_<B>& operator>>=(bitset_base_<B>& b, std::size_t n) noexcept {
         if (n == b.size( )) {
            b.reset( );
         } else if (n != 0) {
            const auto wshift = n / b.bits_per_block( );
            const auto offset = n % b.bits_per_block( );
            if (offset == 0) {
               for (std::size_t i = 0; i <= b.blocks( ) - wshift - 1; ++i) {
                  b.block_begin( )[i] = b.block_begin( )[i + wshift];
               }
            } else {
               const std::size_t sub_offset = b.bits_per_block( ) - offset;
               for (std::size_t i = 0; i < b.blocks( ) - wshift - 1; ++i) {
                  b.block_begin( )[i] = (b.block_begin( )[i + wshift] >> offset) | (b.block_begin( )[i + wshift + 1] << sub_offset);
               }
               b.block_begin( )[b.blocks( ) - wshift - 1] = b.block_begin( )[b.blocks( ) - 1] >> offset;
            }
            std::fill(b.block_begin( ) + b.blocks( ) - wshift, b.block_begin( ) + b.blocks( ), typename bitset_base_<B>::block_type(0));
         }
         return b;
      }

      template<typename B>
      constexpr bitset_base_<B>& operator&=(bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         return std::transform(a.block_begin( ), a.block_end( ), b.block_begin( ), a.block_begin( ), std::bit_and{ }), a;
      }

      template<typename B>
      constexpr bitset_base_<B>& operator|=(bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         return std::transform(a.block_begin( ), a.block_end( ), b.block_begin( ), a.block_begin( ), std::bit_or{ }), a;
      }

      template<typename B>
      constexpr bitset_base_<B>& operator^=(bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         return std::transform(a.block_begin( ), a.block_end( ), b.block_begin( ), a.block_begin( ), std::bit_xor{ }), a;
      }

      template<typename B>
      constexpr bitset_base_<B> operator<<(const bitset_base_<B>& b, std::size_t n) noexcept {
         auto temp = b;
         temp <<= n;
         return temp;
      }

      template<typename B>
      constexpr bitset_base_<B> operator>>(const bitset_base_<B>& b, std::size_t n) noexcept {
         auto temp = b;
         temp >>= n;
         return temp;
      }

      template<typename B>
      constexpr bitset_base_<B> operator&(const bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         auto temp = a;
         temp &= b;
         return temp;
      }

      template<typename B>
      constexpr bitset_base_<B> operator|(const bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         auto temp = a;
         temp |= b;
         return temp;
      }

      template<typename B>
      constexpr bitset_base_<B> operator^(const bitset_base_<B> a, const bitset_base_<B>& b) noexcept {
         auto temp = a;
         temp ^= b;
         return temp;
      }

      template<typename B>
      constexpr bitset_base_<B> operator~(const bitset_base_<B>& a) noexcept {
         auto temp = a;
         temp.flip( );
         return temp;
      }

      template<typename B>
      constexpr bool operator==(const bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         return std::equal(a.block_begin( ), a.block_end( ), b.block_begin( ));
      }

      template<typename B>
      constexpr bool operator!=(const bitset_base_<B>& a, const bitset_base_<B>& b) noexcept {
         return !(a == b);
      }


      template<std::size_t N, typename T = unsigned>
      class static_bitset_ {
      public:
         using block_type = T;

         constexpr static_bitset_( ) noexcept {
            std::for_each(block_data( ), block_data( ) + blocks( ), FUNCTOR(reset_bits));
         }

         static constexpr std::size_t bits_per_block( ) noexcept {
            return bit_size<block_type>( );
         }

         static constexpr std::size_t size( ) noexcept {
            return N;
         }

         static constexpr std::size_t blocks( ) noexcept {
            return ceil_division(size( ), bits_per_block( ));
         }

         static constexpr bool has_trailing_bits( ) noexcept {
            return !is_multiple(size( ), bits_per_block( ));
         }

         static constexpr std::size_t trailing_bits( ) noexcept {
            return multiple_underflow(size( ), bits_per_block( ));
         }

         constexpr T* block_data( ) noexcept {
            return mem_.data( );
         }

         constexpr const T* block_data( ) const noexcept {
            return mem_.data( );
         }

      private:
         std::array<T, ceil_division(N, bit_size<T>( ))> mem_;
      };

      template<typename T = unsigned>
      class adapted_bitset_ {
      public:
         using block_type = T;

         constexpr adapted_bitset_( ) noexcept
         : adapted_bitset_(nullptr, nullptr) {
         }

         constexpr adapted_bitset_(T* i, T* f) noexcept
         : ini_(i), fin_(f) {
         }

         static constexpr std::size_t bits_per_block( ) noexcept {
            return bit_size<block_type>( );
         }

         constexpr std::size_t size( ) const noexcept {
            return (fin_ - ini_) * bit_size<T>( );
         }

         constexpr std::size_t blocks( ) const noexcept {
            return fin_ - ini_;
         }

         constexpr bool has_trailing_bits( ) const noexcept {
            return false;
         }

         constexpr std::size_t trailing_bits( ) const noexcept {
            return 0;
         }

         constexpr T* block_data( ) noexcept {
            return ini_;
         }

         constexpr const T* block_data( ) const noexcept {
            return ini_;
         }

      private:
         T* ini_;
         T* fin_;
      };

      template<typename T = unsigned>
      class const_adapted_bitset_ {
      public:
         using block_type = T;

         constexpr const_adapted_bitset_( ) noexcept
         : const_adapted_bitset_(nullptr, nullptr) {
         }

         constexpr const_adapted_bitset_(const T* i, const T* f) noexcept
         : ini_(i), fin_(f) {
         }

         static constexpr std::size_t bits_per_block( ) noexcept {
            return bit_size<block_type>( );
         }

         constexpr std::size_t size( ) const noexcept {
            return (fin_ - ini_) * bit_size<T>( );
         }

         constexpr std::size_t blocks( ) const noexcept {
            return fin_ - ini_;
         }

         constexpr bool has_trailing_bits( ) const noexcept {
            return false;
         }

         constexpr std::size_t trailing_bits( ) const noexcept {
            return 0;
         }

         constexpr const T* block_data( ) const noexcept {
            return ini_;
         }

      private:
         const T* ini_;
         const T* fin_;
      };

      template<typename T = unsigned, typename A = std::allocator<T>>
      class dynamic_bitset_ : private A {
      public:
         using block_type = T;
         using allocator_type = A;

         constexpr dynamic_bitset_( ) noexcept
         : dynamic_bitset_(A( )) {
         }

         explicit constexpr dynamic_bitset_(A alloc) noexcept
         : A(alloc), mem_(nullptr) {
         }

         explicit constexpr dynamic_bitset_(std::size_t n) noexcept
         : dynamic_bitset_(n, A( )) {
         }

         template<typename RI> //requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
         explicit constexpr dynamic_bitset_(RI ini, RI fin) noexcept
         : dynamic_bitset_(ini, fin, A( )) {
         }

         explicit constexpr dynamic_bitset_(std::size_t n, A alloc) noexcept
         : A(alloc), mem_((n != 0 ? dynamic_construct<detalle_>(alloc, asignar_(n)) : nullptr)) {
            if (mem_ != nullptr) {
               mem_->tam = n;
               std::uninitialized_value_construct_n(block_data( ), blocks( ));
            }
         }

         template<typename RI> //requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
         explicit constexpr dynamic_bitset_(RI ini, RI fin, A alloc) noexcept
         : A(alloc), mem_((ini != fin ? dynamic_construct<detalle_>(alloc, asignar_((fin - ini) * base::bit_size<T>( ))) : nullptr)) {
            if (mem_ != nullptr) {
               mem_->tam = (fin - ini) * base::bit_size<T>( );
               std::uninitialized_copy(ini, fin, block_data( ));
            }
         }

         constexpr dynamic_bitset_(dynamic_bitset_&& d) noexcept
         : A(d.get_allocator( )), mem_(d.mem_) {
            d.mem_ = nullptr;
         }

         constexpr dynamic_bitset_(const dynamic_bitset_& d) noexcept
         : A(d.get_allocator( )), mem_((d.size( ) != 0 ? dynamic_construct<detalle_>(d.get_allocator( ), asignar_(d.size( ))) : nullptr)) {
            if (mem_ != nullptr) {
               mem_->tam = d.size( );
               std::uninitialized_copy_n(d.block_data( ), d.blocks( ), block_data( ));
            }
         }

         ~dynamic_bitset_( ) noexcept {
            libera_( );
         }

         constexpr dynamic_bitset_& operator=(dynamic_bitset_&& d) noexcept {
            return move(*this, d), *this;
         }

         constexpr dynamic_bitset_& operator=(const dynamic_bitset_& d) noexcept {
            return assign(*this, d), *this;
         }

         static constexpr std::size_t bits_per_block( ) noexcept {
            return bit_size<block_type>( );
         }

         constexpr void release( ) noexcept {
            libera_( ), mem_ = nullptr;
         }

         constexpr std::size_t size( ) const noexcept {
            return (mem_ != nullptr ? mem_->tam : 0);
         }

         constexpr std::size_t blocks( ) const noexcept {
            return ceil_division(size( ), bits_per_block( ));
         }

         constexpr bool has_trailing_bits( ) const noexcept {
            return !is_multiple(size( ), bits_per_block( ));
         }

         constexpr std::size_t trailing_bits( ) const noexcept {
            return multiple_underflow(size( ), bits_per_block( ));
         }

         constexpr T* block_data( ) noexcept {
            return (mem_ != nullptr ? mem_->arr : nullptr);
         }

         constexpr const T* block_data( ) const noexcept {
            return (mem_ != nullptr ? mem_->arr : nullptr);
         }

         constexpr auto get_allocator( ) const noexcept {
            return static_cast<const A&>(*this);
         }

      private:
         struct detalle_ {
            std::size_t tam;
            T arr[0];
         };

         constexpr static std::size_t asignar_(std::size_t n) noexcept {
            return ceil_division(sizeof(detalle_) + (ceil_division(n, bits_per_block( )) * sizeof(T)), sizeof(detalle_));
         }

         constexpr void libera_( ) noexcept {
            if (mem_ != nullptr) {
               std::destroy_n(block_data( ), blocks( ));
               dynamic_destroy(get_allocator( ), mem_, asignar_(mem_->tam));
            }
         }

         detalle_* mem_;
      };
   }

   template<std::size_t N, typename T = unsigned>
   using bitset = impl::bitset_base_<impl::static_bitset_<N, T>>;

   template<typename T = unsigned>
   using adapted_bitset = impl::bitset_base_<impl::adapted_bitset_<T>>;

   template<typename T = unsigned>
   using const_adapted_bitset = impl::bitset_base_<impl::const_adapted_bitset_<T>>;

   template<typename T = unsigned, typename A = std::allocator<T>>
   using dynamic_bitset = impl::bitset_base_<impl::dynamic_bitset_<T, A>>;
}

#endif
