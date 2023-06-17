/* zm_export.h */
/* This file was automatically generated. */
/* 2023年  5月 22日 月曜日 08:00:03 JST by zhidao */
#ifndef __ZM_EXPORT_H__
#define __ZM_EXPORT_H__
#include <zeda/zeda_compat.h>
#if defined(__WINDOWS__) && !defined(__CYGWIN__)
# if defined(__ZM_BUILD_DLL__)
#  define __ZM_EXPORT extern __declspec(dllexport)
#  define __ZM_CLASS_EXPORT  __declspec(dllexport)
# else
#  define __ZM_EXPORT extern __declspec(dllimport)
#  define __ZM_CLASS_EXPORT  __declspec(dllimport)
# endif
#else
# define __ZM_EXPORT __EXPORT
# define __ZM_CLASS_EXPORT
#endif
#endif /* __ZM_EXPORT_H__ */
