/*
 * windows.h
 *
 *  Created on: Jan 30, 2019
 *      Author: cgbri
 */

#ifndef __CHEM_NEW_LINK_DLLS_H_
#define __CHEM_NEW_LINK_DLLS_H_

#ifdef __WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#ifdef __cplusplus
namespace windows_link {
#endif

#ifdef __WIN32
typedef HINSTANCE dll_handle_t;
#else
typedef void *dll_handle_t;
#endif

extern dll_handle_t *handles;
extern size_t *num_handles;

/*
 * Call this in your main program to link all the dlls. Works under Windows so far.
 * Pass in the number of libs and an array of null-terminated strings representing the
 * names of the libraries. Will search on the standard path for these, or can use the
 * full path names.
 */
int link_dlls(int num_libs, const char **library_names);

int unlink_dlls(void);

#ifdef __cplusplus
}
#endif

#endif /* WINDOWS_H_ */
