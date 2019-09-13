#ifndef BASE_MEMORY_H
#define BASE_MEMORY_H

#include "bits.h"
#include "integer.h"
#include "intrinsic.h"
#include <algorithm>
#include <climits>
#include <cstddef>
#include <iterator>
#include <memory>
#include <utility>
#include <tuple>
#include <type_traits>

namespace base {
   #ifndef CACHELINE_BYTE_SIZE
      #define CACHELINE_BYTE_SIZE 64
   #endif

   #ifndef L1_CACHE_BYTE_SIZE
      #define L1_CACHE_BYTE_SIZE 32768
   #endif

   #ifndef L2_CACHE_BYTE_SIZE
      #define L2_CACHE_BYTE_SIZE 262144
   #endif

   #ifndef L3_CACHE_BYTE_SIZE
      #define L3_CACHE_BYTE_SIZE 15728640
   #endif

   constexpr std::size_t cacheline_byte_size = CACHELINE_BYTE_SIZE, cacheline_bit_size = cacheline_byte_size * CHAR_BIT;
   constexpr std::size_t l1_cache_byte_size  = L1_CACHE_BYTE_SIZE;
   constexpr std::size_t l2_cache_byte_size  = L2_CACHE_BYTE_SIZE;
   constexpr std::size_t l3_cache_byte_size  = L3_CACHE_BYTE_SIZE;


   template<typename A>
   constexpr bool ignores_deallocation_size_v = false;

   template<typename T>
   constexpr bool ignores_deallocation_size_v<std::allocator<T>> = true;


   template<typename T, typename U>
   constexpr std::size_t member_offset(U T::*m) noexcept {
      return (char*)&(((T*)nullptr)->*m) - (char*)nullptr;
   }


   template<std::size_t N, std::size_t A = floor_pow2(N)>
   struct aligned_storage {
      alignas(A) unsigned char mem[N];
   };

   template<typename T>
   class underlying_type {
   public:
      using value_type = T;

      constexpr value_type& get( ) noexcept {
         return reinterpret_cast<T&>(mem_);
      }

      constexpr const value_type& get( ) const noexcept {
         return reinterpret_cast<const T&>(mem_);
      }

   private:
      aligned_storage<sizeof(T), alignof(T)> mem_;
   };

   template<typename T, std::size_t N> requires alignof(T) >= base::pow2(N)
   class disguised_ptr {
   public:
      constexpr disguised_ptr( ) noexcept
      : disguised_ptr(nullptr) {
      }

      constexpr disguised_ptr(T* p) noexcept
      : p_(reinterpret_cast<std::uintptr_t>(p)) {
      }

      constexpr disguised_ptr(T* p, std::uintptr_t b) noexcept
      : p_(reinterpret_cast<std::uintptr_t>(p) | b) {
      }

      constexpr void disguised_flag(std::uintptr_t b) noexcept {
         *this = disguised_ptr(undisguise( ), b);
      }

      constexpr std::uintptr_t disguised_flag( ) const noexcept {
         return bitmask_until(p_, N);
      }

      constexpr T* undisguise( ) const noexcept {
         return reinterpret_cast<T*>(bitmask_from(p_, N));
      }

   private:
      std::uintptr_t p_;
   };


   template<typename T>
   constexpr auto make_raw_byte_iterator(T& v) noexcept {
      return reinterpret_cast<std::byte*>(&v);
   }

   template<typename T>
   constexpr auto make_raw_byte_iterator(const T& v) noexcept {
      return reinterpret_cast<const std::byte*>(&v);
   }

   template<typename T>
   constexpr auto make_reverse_raw_byte_iterator(T& v) noexcept {
      return std::reverse_iterator<std::byte*>(reinterpret_cast<std::byte*>(&v) + byte_size(v));
   }

   template<typename T>
   constexpr auto make_reverse_raw_byte_iterator(const T& v) noexcept {
      return std::reverse_iterator<const std::byte*>(reinterpret_cast<const std::byte*>(&v) + byte_size(v));
   }

   template<typename T>
   constexpr auto make_most_significant_byte_iterator(T&& v) noexcept {
      if constexpr(std::endian::native == std::endian::big) {
         return make_raw_byte_iterator(std::forward<T>(v));
      } else {
         return make_reverse_raw_byte_iterator(std::forward<T>(v));
      }
   }

   template<typename T>
   constexpr auto make_least_significant_byte_iterator(T&& v) noexcept {
      if constexpr(std::endian::native == std::endian::big) {
         return make_reverse_raw_byte_iterator(std::forward<T>(v));
      } else {
         return make_raw_byte_iterator(std::forward<T>(v));
      }
   }

   using raw_byte_iterator = decltype(make_raw_byte_iterator(std::declval<std::byte>( )));
   using const_raw_byte_iterator = decltype(make_raw_byte_iterator(std::declval<const std::byte>( )));
   using reverse_raw_byte_iterator = decltype(make_reverse_raw_byte_iterator(std::declval<std::byte>( )));
   using const_reverse_raw_byte_iterator = decltype(make_reverse_raw_byte_iterator(std::declval<const std::byte>( )));

   using most_significant_byte_iterator = decltype(make_most_significant_byte_iterator(std::declval<std::byte>( )));
   using const_most_significant_byte_iterator = decltype(make_most_significant_byte_iterator(std::declval<const std::byte>( )));
   using least_significant_byte_iterator = decltype(make_least_significant_byte_iterator(std::declval<std::byte>( )));
   using const_least_significant_byte_iterator = decltype(make_least_significant_byte_iterator(std::declval<const std::byte>( )));


   template<typename T, typename A = std::allocator<T>>
   class unique_ptr;

   template<typename T, typename A>
   class unique_ptr : private A {
   public:
      using value_type = T;
      using allocator_type = A;

      constexpr unique_ptr( ) noexcept
      : unique_ptr(A( )) {
      }

      constexpr unique_ptr(std::nullptr_t) noexcept
      : unique_ptr(A( )) {
      }

      explicit constexpr unique_ptr(A alloc) noexcept
      : A(alloc), mem_(nullptr) {
      }

      explicit constexpr unique_ptr(std::remove_cv_t<T>&& v) noexcept
      : unique_ptr(std::move(v), A( )) {
      }

      explicit constexpr unique_ptr(std::remove_cv_t<T>&& v, A alloc) noexcept
      : A(alloc), mem_(dynamic_construct<std::remove_cv_t<T>>(alloc, 1)) {
         raw_move(*mem_, v);
      }

      explicit constexpr unique_ptr(const std::remove_cv_t<T>& v) noexcept
      : unique_ptr(v, A( )) {
      }

      explicit constexpr unique_ptr(const std::remove_cv_t<T>& v, A alloc) noexcept
      : A(alloc), mem_(dynamic_construct<std::remove_cv_t<T>>(alloc, 1)) {
         raw_assign(*mem_, v);
      }

      constexpr unique_ptr(unique_ptr&& p) noexcept
      : A(p.get_allocator( )), mem_(p.mem_) {
         p.mem_ = nullptr;
      }

      ~unique_ptr( ) noexcept {
         libera_( );
      }

      constexpr unique_ptr(const unique_ptr&) noexcept = delete;

      constexpr unique_ptr& operator=(unique_ptr&& p) noexcept {
         return move(*this, p), *this;
      }

      constexpr unique_ptr& operator=(const unique_ptr&) = delete;

      constexpr T& operator*( ) const noexcept {
         return *mem_;
      }

      constexpr auto operator->( ) const noexcept {
         return get( );
      }

      constexpr auto get( ) const noexcept {
         return mem_;
      }

      constexpr void release( ) noexcept {
         libera_( ), mem_ = nullptr;
      }

      constexpr explicit operator bool( ) const noexcept {
         return mem_ != nullptr;
      }

      constexpr auto get_allocator( ) const noexcept {
         return static_cast<const A&>(*this);
      }

   private:
      constexpr void libera_( ) noexcept {
         if (mem_ != nullptr) {
            dynamic_destroy(get_allocator( ), mem_, 1);
         }
      }

      std::remove_cv_t<T>* mem_;
   };

   namespace impl {
      template<typename T, bool S, typename A = std::allocator<T>>
      class unique_arr_;

      template<typename T, typename A>
      class unique_arr_<T, false, A> : private A {
      public:
         using value_type = T;
         using allocator_type = A;

         constexpr unique_arr_( ) noexcept
         : unique_arr_(A( )) {
         }

         constexpr unique_arr_(std::nullptr_t) noexcept
         : unique_arr_(A( )) {
         }

         explicit constexpr unique_arr_(A alloc) noexcept
         : A(alloc), mem_(nullptr) {
         }

         explicit constexpr unique_arr_(std::size_t n) noexcept
         : unique_arr_(n, A( )) {
         }

         explicit constexpr unique_arr_(std::size_t n, A alloc) noexcept
         : A(alloc), mem_(dynamic_construct<T>(alloc, n)) {
         }

         constexpr unique_arr_(unique_arr_&& p) noexcept
         : A(p.get_allocator( )), mem_(p.mem_) {
            p.mem_ = nullptr;
         }

         ~unique_arr_( ) noexcept {
            libera_( );
         }

         constexpr unique_arr_& operator=(unique_arr_&& p) noexcept {
            return move(*this, p), *this;
         }

         constexpr T& operator*( ) const noexcept {
            return operator[](0);
         }

         constexpr auto operator->( ) const noexcept {
            return get( );
         }

         constexpr value_type& operator[](std::size_t i) const noexcept {
            return mem_[i];
         }

         constexpr auto get( ) const noexcept {
            return mem_;
         }

         constexpr void release( ) noexcept {
            libera_( ), mem_ = nullptr;
         }

         constexpr explicit operator bool( ) const noexcept {
            return mem_ != nullptr;
         }

         constexpr auto get_allocator( ) const noexcept {
            return static_cast<const A&>(*this);
         }

      private:
         constexpr void libera_( ) noexcept {
            if (mem_ != nullptr) {
               dynamic_destroy(get_allocator( ), mem_, 0);
            }
         }

         T* mem_;
      };

      template<typename T, typename A>
      class unique_arr_<T, true, A> : private A {
      public:
         using value_type = T;
         using allocator_type = A;

         constexpr unique_arr_( ) noexcept
         : unique_arr_(A( )) {
         }

         constexpr unique_arr_(std::nullptr_t) noexcept
         : unique_arr_(A( )) {
         }

         explicit constexpr unique_arr_(A alloc) noexcept
         : A(alloc), mem_(nullptr) {
         }

         explicit constexpr unique_arr_(std::size_t n) noexcept
         : unique_arr_(n, A( )) {
         }

         explicit constexpr unique_arr_(std::size_t n, A alloc) noexcept
         : A(alloc), mem_((n != 0 ? dynamic_construct<detalle_>(alloc, asignar_(n)) : nullptr)) {
            if (mem_ != nullptr) {
               std::uninitialized_default_construct_n(mem_->arr, (mem_->tam = n));
            }
         }

         constexpr unique_arr_(unique_arr_&& p) noexcept
         : A(p.get_allocator( )), mem_(p.mem_) {
            p.mem_ = nullptr;
         }

         ~unique_arr_( ) noexcept {
            libera_( );
         }

         constexpr unique_arr_& operator=(unique_arr_&& p) noexcept {
            return move(*this, p), *this;
         }

         constexpr T& operator*( ) const noexcept {
            return operator[](0);
         }

         constexpr auto operator->( ) const noexcept {
            return mem_->arr;
         }

         constexpr value_type& operator[](std::size_t i) const noexcept {
            return mem_->arr[i];
         }

         constexpr T* get( ) const noexcept {
            return (mem_ != nullptr ? mem_->arr : nullptr);
         }

         constexpr void release( ) noexcept {
            libera_( ), mem_ = nullptr;
         }

         constexpr explicit operator bool( ) const noexcept {
            return mem_ != nullptr;
         }

         constexpr auto get_allocator( ) const noexcept {
            return static_cast<const A&>(*this);
         }

      protected:
         constexpr std::size_t size( ) const noexcept {
            return (mem_ != nullptr ? mem_->tam : 0);
         }

      private:
         struct detalle_ {
            std::size_t tam;
            T arr[0];
         };

         constexpr static std::size_t asignar_(std::size_t n) noexcept {
            return ceil_division(sizeof(detalle_) + (n * sizeof(T)), sizeof(detalle_));
         }

         constexpr void libera_( ) noexcept {
            if (mem_ != nullptr) {
               std::destroy_n(mem_->arr, mem_->tam);
               dynamic_destroy(get_allocator( ), mem_, asignar_(mem_->tam));
            }
         }

         detalle_* mem_;
      };
   }

   template<typename T, typename A> requires std::is_trivially_destructible_v<T> && ignores_deallocation_size_v<A>
   struct unique_ptr<T[], A> : public impl::unique_arr_<T, false, A> {
      using impl::unique_arr_<T, false, A>::unique_arr_;
   };

   template<typename T, typename A> requires !(std::is_trivially_destructible_v<T> && ignores_deallocation_size_v<A>)
   struct unique_ptr<T[], A> : public impl::unique_arr_<T, true, A> {
      using impl::unique_arr_<T, true, A>::unique_arr_;
   };

   template<typename T, typename A>
   constexpr bool operator<(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) < b.get( );
   }

   template<typename T, typename A>
   constexpr bool operator<=(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) <= b.get( );
   }

   template<typename T, typename A>
   constexpr bool operator==(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) == b.get( );
   }

   template<typename T, typename A>
   constexpr bool operator>(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) > b.get( );
   }

   template<typename T, typename A>
   constexpr bool operator>=(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) >= b.get( );
   }

   template<typename T, typename A>
   constexpr bool operator!=(const unique_ptr<T, A>& a, const unique_ptr<T, A>& b) {
      return a.get( ) != b.get( );
   }
}

#endif
