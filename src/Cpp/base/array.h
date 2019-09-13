#ifndef BASE_ARRAY_H
#define BASE_ARRAY_H

#include "intrinsic.h"
#include "memory.h"
#include <experimental/array>
#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <utility>
#include <type_traits>

namespace base {
   using std::experimental::make_array;

   namespace impl {
      template<typename T, std::size_t I, std::size_t... S>
      struct matrix {
         using type = std::array<typename matrix<T, S...>::type, I>;
      };

      template<typename T, std::size_t I>
      struct matrix<T, I> {
         using type = std::array<T, I>;
      };
   }

   template<typename T, std::size_t... I>
   using matrix = typename impl::matrix<T, I...>::type;


   template<typename T>
   class array_view {
   public:
      using value_type = T;
      using iterator = const T*;
      using const_iterator = const T*;
      using reverse_iterator = std::reverse_iterator<iterator>;
      using const_reverse_iterator = std::reverse_iterator<const_iterator>;

      constexpr array_view( ) noexcept
      : mem_(nullptr), tam_(0) {
      }

      constexpr array_view(const_iterator ai, const_iterator af) noexcept
      : mem_(ai), tam_(af - ai) {
      }

      template<typename U>
      constexpr array_view(const U& v) noexcept
      : mem_(v.data( )), tam_(v.size( )) {
      }

      template<std::size_t N>
      constexpr array_view(const T(&arr)[N]) noexcept
      : mem_(arr), tam_(N) {
      }

      constexpr array_view(const std::initializer_list<T>& i) noexcept
      : mem_(i.begin( )), tam_(i.size( )) {
      }

      constexpr const T& operator[](std::size_t i) const noexcept {
         return mem_[i];
      }

      constexpr bool empty( ) const noexcept {
         return tam_ == 0;
      }

      constexpr std::size_t size( ) const noexcept {
         return tam_;
      }

      constexpr const T* data( ) const noexcept {
         return mem_;
      }

      constexpr const_iterator begin( ) const noexcept {
         return mem_;
      }

      constexpr const_iterator end( ) const noexcept {
         return mem_ + tam_;
      }

      constexpr const_iterator cbegin( ) const noexcept {
         return begin( );
      }

      constexpr const_iterator cend( ) const noexcept {
         return end( );
      }

      constexpr const_reverse_iterator rbegin( ) const noexcept {
         return const_reverse_iterator(end( ));
      }

      constexpr const_reverse_iterator rend( ) const noexcept {
         return const_reverse_iterator(begin( ));
      }

      constexpr const_reverse_iterator crbegin( ) const noexcept {
         return rbegin( );
      }

      constexpr const_reverse_iterator crend( ) const noexcept {
         return rend( );
      }

      constexpr const value_type& front( ) const noexcept {
         return mem_[0];
      }

      constexpr const value_type& back( ) const noexcept {
         return mem_[tam_ - 1];
      }

   private:
      const T* mem_;
      std::size_t tam_;
   };

   template<typename T>
   constexpr bool operator<(const array_view<T>& a, const array_view<T>& b) noexcept {
      return std::lexicographical_compare(a.begin( ), a.end( ), b.begin( ), b.end( ));
   }

   template<typename T>
   constexpr bool operator==(const array_view<T>& a, const array_view<T>& b) noexcept {
      return (a.size( ) == b.size( ) && std::equal(a.begin( ), a.end( ), b.begin( )));
   }

   template<typename T>
   constexpr bool operator>(const array_view<T>& a, const array_view<T>& b) noexcept {
      return (b < a);
   }

   template<typename T>
   constexpr bool operator<=(const array_view<T>& a, const array_view<T>& b) noexcept {
      return !(a > b);
   }

   template<typename T>
   constexpr bool operator>=(const array_view<T>& a, const array_view<T>& b) noexcept {
      return !(a < b);
   }

   template<typename T>
   constexpr bool operator!=(const array_view<T>& a, const array_view<T>& b) noexcept {
      return !(a == b);
   }


   template<typename T, typename A = std::allocator<std::remove_extent_t<T>>> requires !std::is_array_v<T> && !std::is_const_v<T>
   class dynamic_array: private impl::unique_arr_<T, true, A> {
   public:
      using value_type = typename impl::unique_arr_<T, true, A>::value_type;
      using allocator_type = typename impl::unique_arr_<T, true, A>::allocator_type;
      using iterator = T*;
      using const_iterator = const T*;
      using reverse_iterator = std::reverse_iterator<iterator>;
      using const_reverse_iterator = std::reverse_iterator<const_iterator>;

      using impl::unique_arr_<T, true, A>::unique_arr_;
      using impl::unique_arr_<T, true, A>::size;
      using impl::unique_arr_<T, true, A>::get_allocator;

      template<typename RI>
      constexpr dynamic_array(RI ini, RI fin) noexcept
      : impl::unique_arr_<T, true, A>(fin - ini) {
         std::copy(ini, fin, begin( ));
      }

      constexpr dynamic_array(dynamic_array&&) noexcept = default;
      constexpr dynamic_array(const dynamic_array& p) noexcept
      : impl::unique_arr_<T, true, A>(p.size( ), p.get_allocator( )) {
         std::copy(p.begin( ), p.end( ), begin( ));
      }

      constexpr dynamic_array& operator=(dynamic_array&&) noexcept = default;
      constexpr dynamic_array& operator=(const dynamic_array& p) noexcept {
         return assign(*this, p), *this;
      }

      constexpr bool empty( ) const noexcept {
         return size( ) == 0;
      }

      constexpr T* data( ) noexcept {
         return impl::unique_arr_<T, true, A>::get( );
      }

      constexpr const T* data( ) const noexcept {
         return impl::unique_arr_<T, true, A>::get( );
      }

      constexpr T& operator[](std::size_t i) noexcept {
         return impl::unique_arr_<T, true, A>::operator[](i);
      }

      constexpr const T& operator[](std::size_t i) const noexcept {
         return impl::unique_arr_<T, true, A>::operator[](i);
      }

      constexpr void release( ) noexcept {
         return impl::unique_arr_<T, true, A>::release( );
      }

      constexpr iterator begin( ) noexcept {
         return data( );
      }

      constexpr iterator end( ) noexcept {
         return data( ) + size( );
      }

      constexpr const_iterator begin( ) const noexcept {
         return data( );
      }

      constexpr const_iterator end( ) const noexcept {
         return data( ) + size( );
      }

      constexpr const_iterator cbegin( ) const noexcept {
         return begin( );
      }

      constexpr const_iterator cend( ) const noexcept {
         return end( );
      }

      constexpr reverse_iterator rbegin( ) noexcept {
         return reverse_iterator(end( ));
      }

      constexpr reverse_iterator rend( ) noexcept {
         return reverse_iterator(begin( ));
      }

      constexpr const_reverse_iterator rbegin( ) const noexcept {
         return const_reverse_iterator(end( ));
      }

      constexpr const_reverse_iterator rend( ) const noexcept {
         return const_reverse_iterator(begin( ));
      }

      constexpr const_reverse_iterator crbegin( ) const noexcept {
         return rbegin( );
      }

      constexpr const_reverse_iterator crend( ) const noexcept {
         return rend( );
      }

      constexpr value_type& front( ) noexcept {
         return operator[](0);
      }

      constexpr value_type& back( ) noexcept {
         return operator[](size( ) - 1);
      }

      constexpr const value_type& front( ) const noexcept {
         return operator[](0);
      }

      constexpr const value_type& back( ) const noexcept {
         return operator[](size( ) - 1);
      }
   };

   template<typename T, typename A>
   constexpr bool operator<(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return std::lexicographical_compare(a.begin( ), a.end( ), b.begin( ), b.end( ));
   }

   template<typename T, typename A>
   constexpr bool operator==(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return (a.size( ) == b.size( ) && std::equal(a.begin( ), a.end( ), b.begin( )));
   }

   template<typename T, typename A>
   constexpr bool operator>(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return b < a;
   }

   template<typename T, typename A>
   constexpr bool operator<=(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return !(a > b);
   }

   template<typename T, typename A>
   constexpr bool operator>=(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return !(a < b);
   }

   template<typename T, typename A>
   constexpr bool operator!=(const dynamic_array<T, A>& a, const dynamic_array<T, A>& b) noexcept {
      return !(a == b);
   }
}

#endif
