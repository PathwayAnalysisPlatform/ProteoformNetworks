#ifndef BASE_FUNCTIONAL_H
#define BASE_FUNCTIONAL_H

#include <cstddef>
#include <tuple>
#include <utility>

namespace base {
   #define FUNCTOR(f)                                                       \
      [](auto&&... v) noexcept { return f(std::forward<decltype(v)>(v)...); }

   template<typename... T>
   constexpr void nop(const T&...) noexcept {
      return;
   }

   template<typename T>
   constexpr T& lift(T& v) noexcept {
      return v;
   }

   template<typename T>
   constexpr const T& lift(const T& v) noexcept {
      return v;
   }


   namespace impl {
      template<typename T>
      struct identity_ {
         using type = T;
      };
   }

   template<typename T>
   using type_param = typename impl::identity_<T>;

   template<typename T>
   using identity_t = typename impl::identity_<T>::type;

   template<std::size_t I, typename... T>
   using nth_type_t = std::tuple_element_t<I, std::tuple<T...>>;

   struct empty_type { };


   template<typename F>
   class inverse_predicate {
   public:
      constexpr inverse_predicate( ) noexcept = default;
      constexpr inverse_predicate(F f) noexcept
      : func_(std::move(f)) {
      }

      template<typename T>
      constexpr bool operator()(T&& p) const noexcept {
         return !func_(std::forward<T>(p));
      }

      template<typename T1, typename T2>
      constexpr bool operator()(T1&& p1, T2&& p2) const noexcept {
         return func_(std::forward<T2>(p2), std::forward<T1>(p1));
      }

   private:
      F func_;
   };


   template<typename F>
   class exit_guard {
   public:
      constexpr exit_guard(F&& v) noexcept
      : func_(std::move(v)) {
      }

      constexpr exit_guard(const F& v) noexcept
      : func_(v) {
      }

      ~exit_guard( ) noexcept {
         func_( );
      }

   private:
      F func_;
   };
}

#endif
