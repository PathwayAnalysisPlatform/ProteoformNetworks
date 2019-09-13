#ifndef BASE_ALGORITHM_H
#define BASE_ALGORITHM_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iterator>
#include <utility>
#include <type_traits>

namespace base {
   template<typename II, typename T, typename F>
   constexpr II for_each_n(II i, T n, F&& f) noexcept {
      for (; n != 0; --n, ++i) {
         f(*i);
      }
      return i;
   }


   template<typename RI, typename T, typename P>
   constexpr std::pair<RI, bool> binary_find(RI i, RI f, const T& v, P&& pred) noexcept {
      auto res = std::lower_bound(i, f, v, pred);
      return { res, (res != f && !pred(v, *res)) };
   }

   template<typename RI, typename T>
   constexpr std::pair<RI, bool> binary_find(RI i, RI f, const T& v) noexcept {
      return binary_find(i, f, v, std::less{ });
   }


   template<typename RI1, typename RI2, typename OI, typename P>
   constexpr OI lower_bound(RI1 ai, RI1 af, RI2 bi, RI2 bf, OI wi, P&& pred) noexcept {
      if (ai == af || bi == bf) {
         return std::fill_n(wi, af - ai, bf);
      }

      auto am = ai + (af - bi) / 2;
      auto bj = std::lower_bound(bi, bf, *am, pred);
      wi = lower_bound(ai, am, bi, bj, wi, pred);
      *wi++ = bj;
      wi = lower_bound(am + 1, af, bj, bf, wi, pred);

      return wi;
   }

   template<typename RI1, typename RI2, typename OI>
   constexpr OI lower_bound(RI1 ai, RI1 af, RI2 bi, RI2 bf, OI wi) noexcept {
      return lower_bound(ai, af, bi, bf, wi, std::less{ });
   }


   template<typename T, typename P>
   constexpr int compare_3way(const T& a, const T& b, P&& pred) noexcept {
      return pred(b, a) - pred(a, b);
   }

   template<typename T>
   constexpr int compare_3way(const T& a, const T& b) noexcept {
      return compare_3way(a, b, std::less{ });
   }


   template<typename FI1, typename FI2, typename P>
   constexpr bool memory_compare(FI1 ai, FI1 af, FI2 bi, P&& pred) noexcept {
      for (; ai != af; ++ai, ++bi) {
         if (pred(*ai, *bi)) {
            return true;
         }
         if (pred(*bi, *ai)) {
            return false;
         }
      }
      return false;
   }

   template<typename FI1, typename FI2>
   constexpr bool memory_compare(FI1 ai, FI1 af, FI2 bi) noexcept {
      return memory_compare(ai, af, bi, std::less{ });
   }

   template<typename FI1, typename FI2, typename P>
   constexpr int memory_compare_3way(FI1 ai, FI1 af, FI2 bi, P&& pred) noexcept {
      for (; ai != af; ++ai, ++bi) {
         auto res = compare_3way(*ai, *bi, pred);
         if (res != 0) {
            return res;
         }
      }
      return 0;
   }

   template<typename FI1, typename FI2>
   constexpr int memory_compare_3way(FI1 ai, FI1 af, FI2 bi) noexcept {
      return memory_compare_3way(ai, af, bi, std::less{ });
   }

   template<typename FI1, typename FI2, typename P>
   constexpr int lexicographical_compare_3way(FI1 ai, FI1 af, FI2 bi, FI2 bf, P&& pred) noexcept {
      for (; ai != af; ++ai, ++bi) {
         auto res = compare_3way(*ai, *bi, pred);
         if (res != 0) {
            return res;
         }
      }
      return compare_3way(ai != af, bi != bf);
   }

   template<typename FI1, typename FI2>
   constexpr int lexicographical_compare_3way(FI1 ai, FI1 af, FI2 bi, FI2 bf) noexcept {
      return lexicographical_compare_3way(ai, af, bi, bf, std::less{ });
   }


   namespace impl {
      template<typename FI, typename P>
      constexpr FI find_last_if_(FI ai, FI af, P&& pred, std::forward_iterator_tag) noexcept {
         for (auto res = af; ; ++ai) {
            ai = std::find_if(ai, af, pred);
            if (ai == af) {
               return res;
            }
            res = ai;
         }
      }

      template<typename BI, typename P>
      constexpr BI find_last_if_(BI ai, BI af, P&& pred, std::bidirectional_iterator_tag) noexcept {
         auto ini = std::reverse_iterator(af);
         auto fin = std::reverse_iterator(ai);
         auto res = std::find_if(ini, fin, pred);
         return (res == fin ? af : std::prev(res.base( )));
      }
   }

   template<typename FI, typename P>
   constexpr FI find_last_if(FI ai, FI af, P&& pred) noexcept {
      return impl::find_last_if_(ai, af, pred, typename std::iterator_traits<FI>::iterator_category( ));
   }

   template<typename FI, typename T>
   constexpr FI find_last(FI ai, FI af, const T& v) noexcept {
      return impl::find_last_if_(ai, af, [&](const auto& actual) {
         return actual == v;
      }, typename std::iterator_traits<FI>::iterator_category( ));
   }


   template<typename FI, typename T>
   constexpr FI find_not(FI ai, FI af, const T& v) noexcept {
      return std::find_if_not(ai, af, [&](const auto& actual) {
         return actual == v;
      });
   }


   namespace impl {
      template<typename FI1, typename FI2, typename P>
      constexpr FI1 find_last_of_(FI1 ai, FI1 af, FI2 bi, FI2 bf, P&& pred, std::forward_iterator_tag) noexcept {
         for (auto res = af; ; ++ai) {
            ai = std::find_first_of(ai, af, bi, bf, pred);
            if (ai == af) {
               return res;
            }
            res = ai;
         }
      }

      template<typename BI, typename FI, typename P>
      constexpr BI find_last_of_(BI ai, BI af, FI bi, FI bf, P&& pred, std::bidirectional_iterator_tag) noexcept {
         auto ini = std::reverse_iterator(af);
         auto fin = std::reverse_iterator(ai);
         auto res = std::find_first_of(ini, fin, bi, bf, pred);
         return (res == fin ? af : std::prev(res.base( )));
      }
   }

   template<typename FI1, typename FI2, typename P>
   constexpr FI1 find_last_of(FI1 ai, FI1 af, FI2 bi, FI2 bf, P&& pred) noexcept {
      return impl::find_last_of_(ai, af, bi, bf, pred, typename std::iterator_traits<FI1>::iterator_category( ));
   }

   template<typename FI1, typename FI2>
   constexpr FI1 find_last_of(FI1 ai, FI1 af, FI2 bi, FI2 bf) noexcept {
      return impl::find_last_of_(ai, af, bi, bf, std::equal_to{ }, typename std::iterator_traits<FI1>::iterator_category( ));
   }


   template<typename II, typename OI, typename P>
   constexpr std::pair<II, OI> copy_while(II ai, II af, OI oi, P&& pred) noexcept {
      while (ai != af && pred(*ai)) {
         *oi++ = *ai++;
      }
      return { ai, oi };
   }

   template<typename II, typename OI, typename P>
   constexpr std::pair<II, OI> copy_until(II ai, II af, OI oi, P&& pred) noexcept {
      return copy_while(ai, af, oi, std::not_fn(pred));
   }


   /*namespace impl {
      template<typename FI, typename C> requires C::partition_size( ) > 3
      constexpr std::array<FI, C::partition_size( ) - 1> multiway_partition_(FI ai, FI af, const C& cl, std::forward_iterator_tag) noexcept {
         std::array<std::size_t, C::partition_size( )> cuenta = { };
         for (auto i = ai; i != af; ++i) {
            ++cuenta[cl(*i)];
         }

         auto actual = ai;
         std::array<FI, C::partition_size( )> ini;
         std::array<FI, C::partition_size( ) - 1> fin;
         for (std::size_t i = 0; i < C::partition_size( ) - 1; ++i) {
            ini[i] = actual;
            std::advance(actual, cuenta[i]);
            fin[i] = actual;
         }
         ini.back( ) = actual;

         for (std::size_t i = 0; i < C::partition_size( ) - 1; ++i) {
            for (auto j = ini[i]; j != fin[i]; ) {
               auto offset = cl(*j);
               if (offset == i) {
                  ++j;
               } else {
                  std::swap(*j, *ini[offset]++);
               }
            }
         }

         return fin;
      }

      template<typename FI, typename C> requires C::partition_size( ) == 3
      constexpr std::array<FI, C::partition_size( ) - 1> multiway_partition_(FI ai, FI af, const C& cl, std::forward_iterator_tag) noexcept {
         FI f0 = ai, f1 = ai;
         while (ai != af) {
            auto res = cl(*ai);
            if (res == 0) {
               std::swap(*ai, *f1);
               std::swap(*f0++, *f1++);
            } else if (res == 1) {
               std::swap(*ai, *f1++);
            }
            ++ai;
         }
         return { f0, f1 };
      }

      template<typename BI, typename C> requires C::partition_size( ) == 3
      constexpr std::array<BI, C::partition_size( ) - 1> multiway_partition_(BI ai, BI af, const C& cl, std::bidirectional_iterator_tag) noexcept {
         BI f0 = ai, f1 = af;
         while (ai != f1) {
            auto res = cl(*ai);
            if (res == 2) {
               std::swap(*ai, *--f1);
            } else {
               if (res == 0) {
                  std::swap(*ai, *f0++);
               }
               ++ai;
            }
         }
         return { f0, f1 };
      }

      template<typename FI, typename C, typename IT> requires C::partition_size( ) == 2
      constexpr std::array<FI, C::partition_size( ) - 1> multiway_partition_(FI ai, FI af, const C& cl, IT tag) noexcept {
         return { std::partition(ai, af, std::not_fn(cl)) };
      }

      template<typename FI, typename C, typename IT> requires C::partition_size( ) == 1
      constexpr std::array<FI, C::partition_size( ) - 1> multiway_partition_(FI ai, FI af, const C& cl, IT tag) noexcept {
         return { };
      }
   }

   template<typename FI, typename C>
   constexpr std::array<FI, C::partition_size( ) - 1> multiway_partition(FI ai, FI af, const C& cl) noexcept {
      return impl::multiway_partition_(ai, af, cl, typename std::iterator_traits<FI>::iterator_category( ));
   }*/


   template<typename FI, typename P>
   constexpr FI stable_remove_if(FI ai, FI af, P&& pred) noexcept {
      return std::remove_if(ai, af, pred);
   }

   template<typename FI, typename T>
   constexpr FI stable_remove(FI ai, FI af, const T& v) noexcept {
      return std::remove(ai, af, v);
   }

   namespace impl {
      template<typename BI, typename P>
      constexpr BI remove_if_(BI ai, BI af, P&& pred, std::bidirectional_iterator_tag) noexcept {
         for (;;) {
            auto q = std::find_if(ai, af, pred);
            if (q == af) {
               return q;
            }

            auto s = std::next(q);
            auto m = find_last_if(s, af, std::not_fn(pred), std::bidirectional_iterator_tag( ));
            if (m == af) {
               return q;
            }

            *q = std::move(*m);
            ai = s;
            af = m;
         }
      }

      template<typename FI, typename P>
      constexpr FI remove_if_(FI ai, FI af, P&& pred, std::forward_iterator_tag) noexcept {
         return stable_remove_if(ai, af, pred);
      }
   }

   template<typename FI, typename P>
   constexpr FI remove_if(FI ai, FI af, P&& pred) noexcept {
      return impl::remove_if_(ai, af, pred, typename std::iterator_traits<FI>::iterator_category( ));
   }

   template<typename FI, typename T>
   constexpr FI remove(FI ai, FI af, const T& v) noexcept {
      return impl::remove_if_(ai, af, [&](const auto& actual) {
         return actual == v;
      }, typename std::iterator_traits<FI>::iterator_category( ));
   }


   template<typename T, typename P>
   constexpr void sort2(T& a, T& b, P&& pred) noexcept {
      if (pred(b, a)) {
         std::swap(a, b);
      }
   }

   template<typename T>
   constexpr void sort2(T& a, T& b) noexcept {
      sort2(a, b, std::less{ });
   }


   template<typename RI1, typename RI2>
   constexpr void permute_from(RI1 ai, RI1 af, RI2 bi) noexcept {
      typename std::iterator_traits<RI1>::value_type temp[af - ai];
      for (std::size_t i = 0; i < af - ai; ++i) {
         temp[i] = std::move(ai[bi[i]]);
      }
      std::move(temp, temp + (af - ai), ai);
   }

   template<typename RI1, typename RI2>
   constexpr void permute_to(RI1 ai, RI1 af, RI2 bi) noexcept {
      typename std::iterator_traits<RI1>::value_type temp[af - ai];
      for (std::size_t i = 0; i < af - ai; ++i) {
         temp[bi[i]] = std::move(ai[i]);
      }
      std::move(temp, temp + (af - ai), ai);
   }
}

#endif
