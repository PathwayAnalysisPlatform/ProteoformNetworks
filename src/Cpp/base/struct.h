#ifndef STRUCT_H
#define STRUCT_H

#include "intrinsic.h"
#include "integer.h"
#include "memory.h"
#include <algorithm>
#include <array>
#include <memory>
#include <utility>
#include <tuple>

namespace base {
   template<typename... T>
   struct flexible { };

   template<typename S, typename F, typename A = std::allocator<void>>
   class dynamic_struct;

   template<typename S, typename... T, typename A>
   class dynamic_struct<S, flexible<T[]...>, A> : private A {
   public:
      constexpr dynamic_struct( ) noexcept
      : dynamic_struct(A( )) {
      }

      explicit constexpr dynamic_struct(A alloc) noexcept
      : A(alloc), mem_(nullptr) {
      }

      explicit constexpr dynamic_struct(std::array<std::size_t, sizeof...(T)> tam) noexcept
      : dynamic_struct(S( ), std::move(tam), A( )) {
      }

      explicit constexpr dynamic_struct(std::array<std::size_t, sizeof...(T)> tam, A alloc) noexcept
      : dynamic_struct(S( ), std::move(tam), alloc) {
      }

      explicit constexpr dynamic_struct(S s, std::array<std::size_t, sizeof...(T)> tam) noexcept
      : dynamic_struct(std::move(s), std::move(tam), A( )) {
      }

      explicit constexpr dynamic_struct(S s, std::array<std::size_t, sizeof...(T)> tam, A alloc) noexcept
      : A(alloc), raw_(dynamic_construct<memoria_>(alloc, asignar_((tam2pos_fin_(tam), tam)))) {
         construct(header( ), std::move(s)), mem_->pos_fin = tam;
         unrolled_for<0, sizeof...(T)>([&, this](auto i) {       // "internal error" quitando this
            std::uninitialized_default_construct_n(array<i>( ), size<i>( ));
         });
      }

      constexpr dynamic_struct(dynamic_struct&& s) noexcept
      : A(s.get_allocator( )), raw_(s.raw_) {
         s.raw_ = nullptr;
      }

      constexpr dynamic_struct(const dynamic_struct& s) noexcept
      : A(s.get_allocator( )), raw_((s.empty( ) ? nullptr : dynamic_construct<memoria_>(s.get_allocator( ), asignar_(s.mem_->pos_fin)))) {
         if (mem_ != nullptr) {
            construct(header( ), s.header( )), mem_->pos_fin = s.mem_->pos_fin;
            unrolled_for<0, sizeof...(T)>([&, this](auto i) {  // "internal error" quitando this
               std::uninitialized_copy_n(s.array<i>( ), s.size<i>( ), array<i>( ));
            });
         }
      }

      ~dynamic_struct( ) noexcept {
         libera_( );
      }

      constexpr dynamic_struct& operator=(const dynamic_struct& s) noexcept {
         return assign(*this, s), *this;
      }

      constexpr dynamic_struct& operator=(dynamic_struct&& s) noexcept {
         return move(*this, s), *this;
      }

      constexpr bool empty( ) const noexcept {
         return mem_ == nullptr;
      }

      constexpr void release( ) noexcept {
         libera_( ), mem_ = nullptr;
      }

      constexpr S& header( ) noexcept {
         return std::get<0>(static_cast<base_&>(*mem_));
      }

      constexpr const S& header( ) const noexcept {
         return std::get<0>(static_cast<const base_&>(*mem_));
      }

      template<std::size_t I>
      constexpr std::size_t size( ) const noexcept {
         return (mem_->pos_fin[I] - pos_ini_<I>( )) / sizeof(nth_type_t<I, T...>);
      }

      template<std::size_t I>
      constexpr auto* array( ) noexcept {
         return reinterpret_cast<nth_type_t<I, T...>*>(reinterpret_cast<char*>(mem_) + pos_ini_<I>( ));
      }

      template<std::size_t I>
      constexpr const auto* array( ) const noexcept {
         return reinterpret_cast<const nth_type_t<I, T...>*>(reinterpret_cast<const char*>(mem_) + pos_ini_<I>( ));
      }

      constexpr auto get_allocator( ) const noexcept {
         return static_cast<const A&>(*this);
      }

   private:
      using base_ = std::tuple<S>;
      struct detalle_ : base_ {
         std::array<std::size_t, sizeof...(T)> pos_fin;
      };
      using memoria_ = aligned_storage<std::max({ alignof(detalle_), alignof(T)... })>;

      constexpr static void tam2pos_fin_(std::array<std::size_t, sizeof...(T)>& tam) noexcept {
         std::size_t fin_actual = fin_header_( );
         unrolled_for<0, sizeof...(T)>([&](auto i) {
            fin_actual += multiple_underflow(fin_actual, alignof(nth_type_t<i, T...>));
            fin_actual += tam[i] * sizeof(nth_type_t<i, T...>);
            tam[i] = fin_actual;
         });
      }

      constexpr static std::size_t fin_header_( ) noexcept {
         return member_offset(&detalle_::pos_fin) + sizeof(detalle_::pos_fin);
      }

      constexpr static std::size_t asignar_(const std::array<std::size_t, sizeof...(T)>& pos_fin) noexcept {
         return ceil_division((sizeof...(T) == 0 ? fin_header_( ) : pos_fin.back( )), sizeof(memoria_));
      }

      template<std::size_t I>
      constexpr std::size_t pos_ini_( ) const noexcept {
         std::size_t fin_anterior = (I == 0 ? fin_header_( ) : mem_->pos_fin[I - 1]);
         return fin_anterior + multiple_underflow(fin_anterior, alignof(nth_type_t<I, T...>));
      }

      constexpr void libera_( ) noexcept {
         if (mem_ != nullptr) {
            unrolled_for<0, sizeof...(T)>([&, this](auto i) {     // "internal error" quitando this
               std::destroy_n(array<i>( ), size<i>( ));
            });
            destroy(header( ));
            dynamic_destroy(get_allocator( ), raw_, asignar_(mem_->pos_fin));
         }
      }

      union {
         detalle_* mem_;
         memoria_* raw_;
      };
   };
}

#endif
