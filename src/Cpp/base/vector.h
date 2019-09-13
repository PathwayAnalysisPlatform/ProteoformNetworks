#ifndef VECTOR_H
#define VECTOR_H

#include "intrinsic.h"
#include "memory.h"
#include <array>
#include <atomic>
#include <cstddef>
#include <memory>

namespace base {
   namespace impl {
      template<typename T, std::size_t N, typename ST>
      class fixed_vector_base_ {
      public:
         using value_type = T;
         using iterator = T*;
         using const_iterator = const T*;
         using reverse_iterator = std::reverse_iterator<iterator>;
         using const_reverse_iterator = std::reverse_iterator<const_iterator>;

         constexpr fixed_vector_base_( ) noexcept
         : tam_(0) {
         }

         constexpr fixed_vector_base_(std::size_t n) noexcept
         : tam_(n) {
            std::uninitialized_value_construct_n(begin( ), size( ));
         }

         template<typename RI>
         constexpr fixed_vector_base_(RI i, RI f) noexcept
         : tam_(f - i) {
            std::uninitialized_copy(i, f, begin( ));
         }

         constexpr fixed_vector_base_(fixed_vector_base_&& v)
         : tam_(v.tam_) {
            std::uninitialized_move_n(v.begin( ), v.size( ), begin( )), v.clear( );
         }

         constexpr fixed_vector_base_(const fixed_vector_base_& v)
         : tam_(v.tam_) {
            std::uninitialized_copy_n(v.begin( ), v.size( ), begin( ));
         }

         ~fixed_vector_base_( ) {
            std::destroy_n(begin( ), size( ));
         }

         constexpr fixed_vector_base_& operator=(fixed_vector_base_&& v) {
            return move(*this, v), *this;
         }

         constexpr fixed_vector_base_& operator=(const fixed_vector_base_& v) {
            return assign(*this, v), *this;
         }

         constexpr T& operator[](std::size_t i) noexcept {
            return begin( )[i];
         }

         constexpr const T& operator[](std::size_t i) const noexcept {
            return begin( )[i];
         }

         constexpr bool empty( ) const noexcept {
            return tam_ == 0;
         }

         constexpr std::size_t size( ) const noexcept {
            return tam_;
         }

         static constexpr std::size_t capacity( ) noexcept {
            return N;
         }

         constexpr T* data( ) noexcept {
            return &mem_[0].get( );
         }

         constexpr const T* data( ) const noexcept {
            return &mem_[0].get( );
         }

         constexpr iterator begin( ) noexcept {
            return data( );
         }

         constexpr const_iterator begin( ) const noexcept {
            return data( );
         }

         constexpr const_iterator cbegin( ) const noexcept {
            return data( );
         }

         constexpr iterator end( ) noexcept {
            return data( ) + size( );
         }

         constexpr const_iterator end( ) const noexcept {
            return data( ) + size( );
         }

         constexpr const_iterator cend( ) const noexcept {
            return data( ) + size( );
         }

         constexpr reverse_iterator rbegin( ) noexcept {
            return reverse_iterator(end( ));
         }

         constexpr const_reverse_iterator rbegin( ) const noexcept {
            return reverse_iterator(end( ));
         }

         constexpr const_reverse_iterator crbegin( ) const noexcept {
            return reverse_iterator(end( ));
         }

         constexpr reverse_iterator rend( ) noexcept {
            return reverse_iterator(begin( ));
         }

         constexpr const_reverse_iterator rend( ) const noexcept {
            return const_reverse_iterator(begin( ));
         }

         constexpr const_reverse_iterator crend( ) const noexcept {
            return const_reverse_iterator(begin( ));
         }

         constexpr value_type& front( ) noexcept {
            return *begin( );
         }

         constexpr const value_type& front( ) const noexcept {
            return *begin( );
         }

         constexpr value_type& back( ) noexcept {
            return *rbegin( );
         }

         constexpr const value_type& back( ) const noexcept {
            return *rbegin( );
         }

         void push_back(T&& v) {
            construct(data( )[tam_++], std::move(v));
         }

         void push_back(const T& v) {
            construct(data( )[tam_++], v);
         }

         template<typename... P>
         constexpr void emplace_back(P&&... v) {
            construct(data( )[tam_++], std::forward<P>(v)...);
         }

         constexpr void pop_back( ) {
            destroy(data( )[--tam_]);
         }

         constexpr void resize(std::size_t n) {
            if (n < size( )) {
               std::destroy(begin( ) + n, begin( ) + size( ));
            } else if (n > size( )) {
               std::uninitialized_value_construct(begin( ) + size( ), begin( ) + n);
            }
            tam_ = n;
         }

         void clear( ) {
            resize(0);
         }

      private:
         ST tam_;
         std::array<underlying_type<T>, N> mem_;
      };
   }

   template<typename T, std::size_t N>
   class fixed_vector : public impl::fixed_vector_base_<T, N, std::size_t> {
      using impl::fixed_vector_base_<T, N, std::size_t>::fixed_vector_base_;
   };

   template<typename T, std::size_t N>
   class concurrent_fixed_vector : public impl::fixed_vector_base_<T, N, std::atomic<std::size_t>> {
      using impl::fixed_vector_base_<T, N, std::atomic<std::size_t>>::fixed_vector_base_;
      constexpr void pop_back( ) = delete;
   };
}

#endif
