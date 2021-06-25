
#pragma once

#ifdef TLS_EXPORTS
#    define TLS_API __declspec(dllexport)
#else
#    define TLS_API __declspec(dllimport)
#endif