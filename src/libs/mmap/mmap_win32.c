#define WIN32_LEAN_AND_MEAN
#include "windows.h"

struct MMapContext
{
	HANDLE file_handle;
	HANDLE mmap_handle;
	size_t len;
} ;

struct MMapContext *mmap_file_read(const char *filename)
{
	HANDLE file_handle = CreateFileA(filename, GENERIC_READ, FILE_SHARE_READ,
		NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (file_handle == INVALID_HANDLE_VALUE) {
		return 0;
	}
	LARGE_INTEGER file_size = { 0 };
	GetFileSizeEx(file_handle, &file_size);
	HANDLE mmap_handle = CreateFileMappingA(file_handle, NULL, PAGE_READONLY, 0, 0, NULL);
	struct MMapContext *ret = HeapAlloc(GetProcessHeap(), HEAP_GENERATE_EXCEPTIONS, sizeof(struct MMapContext));
	ret->file_handle = file_handle;
	ret->mmap_handle = mmap_handle;
	ret->len = file_size.QuadPart;
	return ret;
}

struct MMapContext *mmap_file_write(const char *filename, size_t file_size)
{
	HANDLE file_handle = CreateFileA(filename, GENERIC_WRITE|GENERIC_READ, FILE_SHARE_READ,
		NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (file_handle == INVALID_HANDLE_VALUE) {
		return 0;
	}
	LARGE_INTEGER s;
	s.QuadPart = file_size;
	HANDLE mmap_handle = CreateFileMappingA(file_handle, NULL, PAGE_READWRITE, s.HighPart, s.LowPart, NULL);
	struct MMapContext *ret = HeapAlloc(GetProcessHeap(), HEAP_GENERATE_EXCEPTIONS, sizeof(struct MMapContext));
	ret->file_handle = file_handle;
	ret->mmap_handle = mmap_handle;
	ret->len = file_size;
	return ret;
}

void mmap_close(struct MMapContext *context)
{
	CloseHandle(context->mmap_handle);
	CloseHandle(context->file_handle);
	HeapFree(GetProcessHeap(), 0, context);
}

size_t mmap_size(struct MMapContext* context)
{
	return context->len;
}

void *mmap_get_ptr_read(struct MMapContext *context)
{
	void * ret = MapViewOfFile(context->mmap_handle, FILE_MAP_READ, 0, 0, 0);
	return ret;
}

void *mmap_get_ptr_write(struct MMapContext *context)
{
	void * ret = MapViewOfFile(context->mmap_handle, FILE_MAP_WRITE, 0, 0, 0);
	return ret;
}
