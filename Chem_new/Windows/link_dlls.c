/*
 * link_dlls.c
 *
 *  Created on: Jan 30, 2019
 *      Author: cgbri
 */

#include "link_dlls.h"

#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <stddef.h>
#include <stdio.h>

#ifdef __WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

dll_handle_t *handles = NULL;
size_t *num_handles = NULL;

int link_dlls(int num_libs, const char **libs) {
	char *is_linked = calloc(num_libs, sizeof(char));
	int num_linked = 0;
	if(handles == NULL) {
		handles = malloc(num_libs * sizeof(dll_handle_t));
		num_handles = malloc(sizeof(int));
		*num_handles = 0;
	} else {
		handles = realloc(handles, *num_handles * sizeof(dll_handle_t));
	}

	for(int i = 0; i < num_libs; i++) {
		if(is_linked[i] == 0) {
#ifdef __WIN32
			dll_handle_t handle = LoadLibrary(libs[i]);
#else
			dll_handle_t handle = dlopen(libs[i], RTLD_NOW);
			if(handle != NULL) {
				is_linked[i] = 1;
				num_linked++;
				handles[i + *num_handles] = handle;
			} else {
#ifdef ELIBACC
				errno = ELIBACC;
#else
				errno = EACCES;
#endif
#ifdef __cplusplus
				throw(std::runtine_exception("Error while loading libraries!"));
#else
				perror("Error while loading libraries!\n");
				*num_handles = num_linked;
				return (-1);
#endif
			}
		}
	}

	*num_handles = num_linked;
	free(is_linked);
	return (num_linked);
}

int unlink_dlls(void) {
	if(num_handles == NULL || handles == NULL) {
		return (-1);
	}
	for(int i = 0; i < *num_handles; i++) {
#ifdef __WIN32
		FreeLibrary(handles[i]);
#else
		dlclose(handles[i]);
#endif
	}
	free(handles);
	free(num_handles);
	handles = NULL;
	num_handles = NULL;
	return (1);
}
