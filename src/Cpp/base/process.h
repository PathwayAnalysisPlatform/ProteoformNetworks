#ifndef BASE_PROCESS_H
#define BASE_PROCESS_H

#include <algorithm>

#ifdef _WIN32
   #define WIN32_LEAN_AND_MEAN
   #include <windows.h>
#else
   #include <limits.h>
   #include <unistd.h>
#endif

namespace base {
   template<typename OI>
   OI process_path(OI salida) noexcept {
   #ifdef _WIN32
      wchar_t buffer[MAX_PATH];
      return std::copy(buffer, buffer + ::GetModuleFileNameW(nullptr, buffer, MAX_PATH), salida);
   #else
      char buffer[PATH_MAX];
      return std::copy(buffer, buffer + std::max(readlink("/proc/self/exe", buffer, PATH_MAX), ssize_t(0)), salida);
   #endif
   }
}

#endif
