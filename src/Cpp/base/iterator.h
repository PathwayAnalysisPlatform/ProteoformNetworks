#ifndef BASE_ITERATOR_H
#define BASE_ITERATOR_H

#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

namespace base {
   template<typename T> //requires std::is_integral_v<T>
   class integer_iterator {
   public:
      using iterator_category = std::random_access_iterator_tag;
      using value_type = T;
      using difference_type = decltype(std::declval<T&>( ) - std::declval<T&>( ));
      using pointer = T*;
      using reference = T&;

      integer_iterator( ) noexcept = default;

      integer_iterator(T v) noexcept
      : v_(std::move(v)) {
      }

      constexpr T& operator*( ) noexcept {
         return v_;
      }

      constexpr const T& operator*( ) const noexcept {
         return v_;
      }

      constexpr integer_iterator& operator++( ) noexcept {
         return ++v_, *this;
      }

      constexpr integer_iterator operator++(int) noexcept {
         auto temp = *this;
         return ++v_, temp;
      }

      constexpr integer_iterator& operator--( ) noexcept {
         return --v_, *this;
      }

      constexpr integer_iterator operator--(int) noexcept {
         auto temp = *this;
         return --v_, temp;
      }

      constexpr integer_iterator& operator+=(const T& d) noexcept {
         return v_ += d, *this;
      }

      constexpr integer_iterator& operator-=(const T& d) noexcept {
         return v_ += d, *this;
      }

   private:
      T v_;
   };

   template<typename T>
   constexpr auto operator-(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return *a - *b;
   }

   template<typename T>
   constexpr bool operator<(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return *a < *b;
   }

   template<typename T>
   constexpr bool operator==(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return *a == *b;
   }

   template<typename T>
   constexpr bool operator>(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return (b < a);
   }

   template<typename T>
   constexpr bool operator<=(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return !(a > b);
   }

   template<typename T>
   constexpr bool operator>=(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return !(a < b);
   }

   template<typename T>
   constexpr bool operator!=(const integer_iterator<T>& a, const integer_iterator<T>& b) noexcept {
      return !(a == b);
   }

   namespace impl {
      template<typename II>
      class range_ {
      public:
         constexpr range_(II i, II f) noexcept
         : i_(i), f_(f) {
         }

         constexpr bool empty( ) const noexcept {
            return i_ == f_;
         }

         constexpr std::size_t size( ) const noexcept {
            return std::distance(i_, f_);
         }

         constexpr II begin( ) const noexcept {
            return i_;
         }

         constexpr II end( ) const noexcept {
            return f_;
         }

         void advance_begin(std::ptrdiff_t n) noexcept {
            std::advance(i_, n);
         }

         void advance_end(std::ptrdiff_t n) noexcept {
            std::advance(f_, n);
         }

         constexpr auto& operator[](std::size_t i) const noexcept {
            return begin( )[i];
         }

      private:
         II i_, f_;
      };
   }

   template<typename II>
   class range : public impl::range_<II> {
   public:
      template<typename T>
      constexpr range(T& v) noexcept
      : impl::range_<II>(std::begin(v), std::end(v)) {
      }

      constexpr range(II i, II f) noexcept
      : impl::range_<II>(i, f) {
      }
   };

   template<typename II>
   class reverse_range : public impl::range_<std::reverse_iterator<II>> {
   public:
      template<typename T>
      constexpr reverse_range(T& v) noexcept
      : impl::range_<std::reverse_iterator<II>>(std::rbegin(v), std::rend(v)) {
      }

      constexpr reverse_range(II i, II f) noexcept
      : impl::range_<std::reverse_iterator<II>>(std::reverse_iterator(f), std::reverse_iterator(i)) {
      }
   };

   template<typename II> //requires !std::is_integral_v<II>
   range(II, II) -> range<II>;

   template<typename T> //requires std::is_integral_v<T>
   range(const T&, const T&) -> range<integer_iterator<T>>;

   template<typename T>
   range(T& v) -> range<decltype(std::begin(v))>;

   template<typename II> //requires !std::is_integral_v<II>
   reverse_range(II, II) -> reverse_range<II>;

   template<typename T> //requires std::is_integral_v<T>
   reverse_range(const T&, const T&) -> reverse_range<integer_iterator<T>>;

   template<typename T>
   reverse_range(T& v) -> reverse_range<decltype(std::begin(v))>;
}

#endif
