#ifndef BASE_COMPRESSED_BITSTREAM_H
#define BASE_COMPRESSED_BITSTREAM_H

#include "bits.h"
#include "bitset.h"
#include "integer.h"
#include "iterator.h"
#include "struct.h"
#include "memory.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>
#include <type_traits>

namespace base {
   namespace impl {
      constexpr std::pair<std::size_t, std::uint8_t> memoria_maxima_(std::size_t capacidad) noexcept {
         std::size_t nodos = 0, bytes = 0;
         if (capacidad != 0) {
            do {
               capacidad = ceil_division(capacidad, 256);
               nodos += capacidad, bytes += 1;
            } while (capacidad != 1);
         }
         return { nodos, bytes };
      }


      template<typename RI>
      constexpr auto crea_bitset_(const std::pair<RI, RI>& rango) noexcept {
         bitset<256, typename std::iterator_traits<RI>::value_type> b;
         std::copy(rango.first, rango.second, b.block_begin( ));
         return b;
      }

      template<typename T>
      constexpr std::size_t bloques_segmento_(std::uint8_t octeto) noexcept {
         return int_pow(256, octeto) * (256 / bit_size<T>( ));
      }

      template<typename RI>
      constexpr std::pair<RI, RI> rango_segmento_(std::uint8_t octeto, std::size_t pos, RI pi, RI pf) noexcept {
         return { pi + (bloques_segmento_<typename std::iterator_traits<RI>::value_type>(octeto) * pos), pi + std::min(bloques_segmento_<typename std::iterator_traits<RI>::value_type>(octeto) * (pos + 1), std::size_t(pf - pi)) };
      }


      enum tipo_ { LISTA = 0, BITSET = 1 } ;
      enum gloton_ { CERO_GLOTON = 0, UNO_GLOTON = 1 } ;
      enum listados_ { CEROS_LISTADOS = 0, UNOS_LISTADOS = 1 } ;
      struct header_ {
         std::uint8_t tipo : 1 = { };
         std::uint8_t gloton : 1 = { };
         std::uint8_t listados : 1 = { };
         std::uint8_t tamanyo : 5 = { };
      };

      template<typename T>
      struct estadistica_nodo_ {
         constexpr estadistica_nodo_( ) noexcept {
            gloton.set( );
         }

         constexpr void registra(std::size_t i, const std::array<unsigned, 2>& v) noexcept {
            visitado[i] = true;
            poblado[i] = (v[0] != 256);
            gloton[i] = (v[0] == 256 || v[1] == 256);
         }

         bitset<256, T> visitado, poblado, gloton;
      };

      struct uso_memoria_ {
         uso_memoria_( ) noexcept
         : bits_headers(0), bytes_listas(0), bitsets(0) {
         }

         std::size_t bits_headers;
         std::size_t bytes_listas, bitsets;
      };

      template<typename T>
      constexpr std::array<unsigned, 2> preprocesa_hoja_(header_& h, const bitset<256, T>& b, uso_memoria_& m) noexcept {
         unsigned existen = b.count( ), faltan = 256 - existen;
         if (existen < 32 || faltan < 32) {
            h = { .tipo = LISTA, .listados = (existen < faltan ? UNOS_LISTADOS : CEROS_LISTADOS), .tamanyo = std::uint8_t(std::min(existen, faltan)) };
            m.bits_headers += 2 + 5, m.bytes_listas += h.tamanyo;
         } else {
            h = { .tipo = BITSET };
            m.bits_headers += 1, m.bitsets += 1;
         }

         return { faltan, existen };
      }

      template<typename T>
      constexpr std::array<unsigned, 2> preprocesa_interno_(header_& h, bitset<256, T>& s, const estadistica_nodo_<T>& e, uso_memoria_& m, bool hoja_abajo) noexcept {
         std::array<unsigned, 2> cuenta_gloton = { 0, 0 };
         e.gloton.visit_set([&](auto i) {
            ++cuenta_gloton[e.poblado[i]];
         });

         if (cuenta_gloton[0] >= cuenta_gloton[1]) {
            preprocesa_hoja_(h, s = e.poblado, m);
            h.gloton = CERO_GLOTON;
            m.bits_headers += 1 - ((e.visitado & e.gloton & ~e.poblado).count( ) * (2 + !hoja_abajo + 5));
         } else {
            preprocesa_hoja_(h, s = e.poblado & e.gloton, m);
            h.gloton = UNO_GLOTON;
            m.bits_headers += 1 - ((e.visitado & e.gloton &  e.poblado).count( ) * (2 + !hoja_abajo + 5));
         }

         return cuenta_gloton;
      }

      template<typename T, typename RI> requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
      constexpr std::array<unsigned, 2> preprocesa_(header_** htorre, bitset<256, T>** storre, std::uint8_t octeto, std::size_t pos, RI pi, RI pf, uso_memoria_& m) noexcept {
         if (octeto == 0) {
            return preprocesa_hoja_(htorre[octeto][pos], crea_bitset_(rango_segmento_(octeto, pos, pi, pf)), m);
         } else {
            estadistica_nodo_<T> e;
            auto rango = rango_segmento_(octeto, pos, pi, pf);
            for (auto i : range(std::size_t(0), ceil_division(rango.second - rango.first, bloques_segmento_<T>(octeto - 1)))) {
               e.registra(i, preprocesa_(htorre, storre, octeto - 1, 256 * pos + i, pi, pf, m));
            }
            return preprocesa_interno_(htorre[octeto][pos], storre[octeto][pos], e, m, octeto == 1);
         }
      }

      template<typename T, typename RI> requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
      constexpr uso_memoria_ preprocesa_(header_** htorre, bitset<256, T>** storre, std::uint8_t bytes, RI pi, RI pf) noexcept {
         uso_memoria_ m;
         preprocesa_(htorre, storre, bytes - 1, 0, pi, pf, m);
         return m;
      }


      static_assert(std::endian::native == std::endian::little);

      constexpr std::uint8_t lee_bits_(std::pair<const std::uint8_t*, std::uint8_t>& datos, std::uint8_t n) noexcept {
         auto res = get_n_bits(*reinterpret_cast<const std::uint16_t*>(datos.first), datos.second, n);
         datos.second += n;
         datos.first += (datos.second >= 8);
         datos.second %= 8;
         return res;
      }

      constexpr void escribe_bits_(std::pair<std::uint8_t*, std::uint8_t>& datos, std::uint8_t n, std::uint8_t v) noexcept {
         write_n_bits(*reinterpret_cast<std::uint16_t*>(datos.first), datos.second, n, v);
         datos.second += n;
         datos.first += (datos.second >= 8);
         datos.second %= 8;
      }


      template<typename T, typename RI> requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
      constexpr void comprime_(header_** htorre, bitset<256, T>** storre, std::uint8_t octeto, std::size_t pos, RI pi, RI pf, std::pair<std::uint8_t*, std::uint8_t>& hw, std::uint8_t*& lw, bitset<256, T>*& bw) noexcept {
         auto h = htorre[octeto][pos];
         auto s = (octeto != 0 ? storre[octeto][pos] : crea_bitset_(rango_segmento_(octeto, pos, pi, pf)));
         escribe_bits_(hw, 1, h.tipo);
         escribe_bits_(hw, bool(octeto), h.gloton);
         if (h.tipo == LISTA) {
            escribe_bits_(hw, 1, h.listados);
            escribe_bits_(hw, 5, h.tamanyo);
            (h.listados == UNOS_LISTADOS ? s : ~s).visit_set([&](auto i) {
               *lw++ = i;
            });
         } else {
            *bw++ = s;
         }

         if (octeto != 0) {
            (h.gloton == CERO_GLOTON ? s : ~s).visit_set([&](auto i) {
               comprime_(htorre, storre, octeto - 1, 256 * pos + i, pi, pf, hw, lw, bw);
            });
         }
      }

      template<typename T, typename RI> requires std::is_same_v<T, typename std::iterator_traits<RI>::value_type>
      constexpr void comprime_(header_** htorre, bitset<256, T>** storre, std::uint8_t bytes, RI pi, RI pf, std::uint8_t* hw, std::uint8_t* lw, bitset<256, T>* bw) noexcept {
         std::pair<std::uint8_t*, std::uint8_t> hwc = { hw, 0 };
         comprime_(htorre, storre, bytes - 1, 0, pi, pf, hwc, lw, bw);
      }


      template<typename T, typename F>
      constexpr void descomprime_(std::pair<const std::uint8_t*, std::uint8_t>& hi, const std::uint8_t*& li, const bitset<256, T>*& bi, std::uint8_t octeto, std::size_t prefix, F&& f) noexcept {
         auto h = header_{ };
         auto s = bitset<256, T>( );
         h.tipo = lee_bits_(hi, 1);
         h.gloton = lee_bits_(hi, bool(octeto));
         if (h.tipo == LISTA) {
            h.listados = lee_bits_(hi, 1);
            h.tamanyo = lee_bits_(hi, 5);
            li = for_each_n(li, h.tamanyo, [&](auto i) {
               s[i] = true;
            }), s = (h.listados == UNOS_LISTADOS ? s : ~s);
         } else {
            s = *bi++;
         }

         if (auto q = 0; octeto != 0) {
            auto no_visitado = [&](auto ini, auto fin) {
               if (ini != fin && h.gloton == UNO_GLOTON) {
                  f(prefix | (ini << (8 * octeto)), prefix | (fin << (8 * octeto)));
               }
            };
            (h.gloton == CERO_GLOTON ? s : ~s).visit_set([&](auto i) {
               no_visitado(q, i), q = i + 1;
               descomprime_(hi, li, bi, octeto - 1, prefix | (i << (8 * octeto)), f);
            }), no_visitado(q, 256);
         } else {
            f(prefix, s);
         }
      }

      template<typename T, typename F>
      constexpr void descomprime_(const std::uint8_t* hi, const std::uint8_t* li, const bitset<256, T>* bi, std::uint8_t bytes, F&& f) noexcept {
         std::pair<const std::uint8_t*, std::uint8_t> hic = { hi, 0 };
         descomprime_(hic, li, bi, bytes - 1, 0, f);
      }


      template<typename F>
      class visitante_set_ {
      public:
         constexpr visitante_set_(F& f) noexcept
         : f_(f) {
         }

         constexpr void operator()(std::size_t ini, std::size_t fin) noexcept {
            for (auto i = ini; i != fin; ++i) {
               f_(i);
            }
         }

         template<typename T>
         constexpr void operator()(std::size_t prefix, const bitset<256, T>& b) noexcept {
            b.visit_set([&](auto i) {
               f_(prefix + i);
            });
         }

      private:
         F& f_;
      };
   }

   template<typename T = unsigned, typename A = std::allocator<void>> requires is_pow2(bit_size<T>( )) && bit_size<T>( ) <= 256
   class compressed_bitstream {
   public:
      using block_type = T;

      constexpr compressed_bitstream( ) noexcept = default;

      template<typename RI> requires std::is_same_v<block_type, typename std::iterator_traits<RI>::value_type>
      constexpr compressed_bitstream(RI pi, RI pf) noexcept
      : compressed_bitstream(pi, pf, A( )) {
      }

      template<typename RI> requires std::is_same_v<block_type, typename std::iterator_traits<RI>::value_type>
      constexpr compressed_bitstream(RI pi, RI pf, A alloc) noexcept {
         if (adapted_bitset<const block_type>(pi, pf).any( )) {
            auto capacidad = adapted_bitset<const block_type>(pi, pf).size( );
            auto [max_nodos, bytes] = impl::memoria_maxima_(capacidad);

            impl::header_ buffer_headers[max_nodos], *htorre[bytes], *hactual = buffer_headers;
            bitset<256, block_type> buffer_segmentos[max_nodos - ceil_division(capacidad, 256)], *storre[bytes], *sactual = buffer_segmentos;
            for (std::uint8_t i = 0; i < bytes; ++i) {
               htorre[i] = hactual, hactual += ceil_division(capacidad, int_pow(std::size_t(256), i + 1));
               storre[i] = sactual, sactual += ceil_division(capacidad, int_pow(std::size_t(256), i + 1)) * (i != 0);
            }

            impl::uso_memoria_ m = impl::preprocesa_(htorre, storre, bytes, pi, pf);
            s_ = decltype(s_)(bytes, { ceil_division(m.bits_headers, 8), m.bytes_listas, m.bitsets });
            std::fill_n(s_.template array<0>( ), s_.template size<0>( ), 0);
            impl::comprime_(htorre, storre, bytes, pi, pf, s_.template array<0>( ), s_.template array<1>( ), s_.template array<2>( ));
         }
      }

      constexpr bool none( ) const noexcept {
         return s_.empty( );
      }

      constexpr bool has_resource( ) const noexcept {
         return !s_.empty( );
      }

      constexpr void release( ) noexcept {
         s_.release( );
      }

      template<typename F> requires std::is_invocable_v<F, std::size_t>
      constexpr void visit_set(F&& f) const noexcept {
         if (!none( )) {
            impl::descomprime_(s_.template array<0>( ), s_.template array<1>( ), s_.template array<2>( ), s_.header( ), impl::visitante_set_(f));
         }
      }

      template<typename F> requires std::is_invocable_v<F, std::size_t, std::size_t> && std::is_invocable_v<F, std::size_t, const bitset<256, block_type>&>
      constexpr void visit_set(F&& f) const noexcept {
         if (!none( )) {
            impl::descomprime_(s_.template array<0>( ), s_.template array<1>( ), s_.template array<2>( ), s_.header( ), f);
         }
      }

   private:
      dynamic_struct<std::uint8_t, flexible<std::uint8_t[], std::uint8_t[], bitset<256, block_type>[]>, A> s_;
   };
}

#endif
