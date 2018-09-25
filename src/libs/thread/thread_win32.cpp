#define WIN32_LEAN_AND_MEAN 
#include <windows.h>

struct ThreadHandle{
    HANDLE handle;
};
struct Semaphore{
    HANDLE handle;
};

extern "C"
struct ThreadHandle *thread_start(unsigned long (*func)(void *),void *param)
{
    struct ThreadHandle *ret = new ThreadHandle;
    ret->handle = CreateThread(NULL,0,func,param,0,NULL);
    return ret;
}

extern "C"
void thread_wait(struct ThreadHandle *handle)
{
    WaitForSingleObject(handle->handle,INFINITE);
    delete handle;
}

extern "C"
struct Semaphore *thread_semaphore_create(unsigned int val, const char *name)
{
    struct Semaphore *ret = new Semaphore;
    ret->handle =CreateSemaphoreA(0,val,100000,name);
    return ret;
}

extern "C"
void thread_semaphore_destroy(struct Semaphore *sem)
{
    CloseHandle(sem->handle);
    delete sem;
}

extern "C"
void thread_semaphore_post(struct Semaphore *sem)
{
    ReleaseSemaphore(sem->handle,1,0);
}

extern "C"
void thread_semaphore_wait(struct Semaphore *sem)
{
    WaitForSingleObject(sem->handle,INFINITE);
}
