#define WIN32_LEAN_AND_MEAN
#include "windows.h"

struct MMapContext
{
	HANDLE file_handle;
	HANDLE mmap_handle;
} ;

struct MMapContext *mmap_file_read(const char *filename)
{
	HANDLE file_handle = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ,
		NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	HANDLE mmap_handle = CreateFileMappingA(file_handle, NULL, PAGE_READONLY, 0, 0, NULL);
	struct MMapContext *ret = HeapAlloc(GetProcessHeap(), HEAP_GENERATE_EXCEPTIONS, sizeof(struct MMapContext));
	ret->file_handle = file_handle;
	ret->mmap_handle = mmap_handle;
	return ret;
}

void mmap_close(struct MMapContext *context)
{
	CloseHandle(context->mmap_handle);
	CloseHandle(context->file_handle);
	HeapFree(GetProcessHeap(), 0, context);
}

void *mmap_get_ptr(struct MMapContext *context)
{
	void * ret = MapViewOfFile(context->mmap_handle, FILE_MAP_READ, 0, 0, 0);
	return ret;
}