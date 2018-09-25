#include "thread_pthread.h"
#include "thread.h"
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>

typedef struct
{
    unsigned long (*func)(void *);
    void *param;
} ThreadWrapperData;

static void * thread_wrapper(void *data){
    ThreadWrapperData *thread_data = (ThreadWrapperData*)data;
    thread_data->func(thread_data->param);
    free(thread_data);
    return NULL;
}

struct ThreadHandle *thread_start(unsigned long (*func)(void *),void *param)
{
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr,8388608);
    ThreadWrapperData *thread_data =
        (ThreadWrapperData *)malloc(sizeof(ThreadWrapperData));
    thread_data->func = func;
    thread_data->param = param;
    pthread_t handle;
    pthread_create(&handle,&attr,thread_wrapper,(void*)thread_data);
    struct ThreadHandle *ret = calloc(1,sizeof(struct ThreadHandle));
    ret->handle = handle;
    return ret;
}

void thread_wait(struct ThreadHandle *handle)
{
    pthread_join(handle->handle,NULL);
}

//TODO(Vidar):sem_init is depricated on osx, use sem_open instead...
struct Semaphore *thread_semaphore_create(unsigned int val, const char *name)
{
    struct Semaphore *ret = calloc(1,sizeof(struct Semaphore));
    ret->semaphore = sem_open(name, O_CREAT | O_EXCL, 0644, val);
    if(ret->semaphore == SEM_FAILED){
        printf("Error opening semaphore: %d\n", errno);
        switch(errno){
            case EACCES: printf(" The semaphore exists, but the caller does not have permission to open it.\n");
                         break;

            case EEXIST: printf(" Both O_CREAT and O_EXCL were specified in oflag, but a semaphore with this name already exists.\n");
                         break;

            case EINVAL: printf(" value was greater than SEM_VALUE_MAX. or name consists of just \"/\", followed by no other characters.  \n");
                         break;
            case EMFILE: printf(" The per-process limit on the number of open file descriptors has been reached.\n");
                         break;

            case ENAMETOOLONG: printf("name was too long.\n");
                         break;

            case ENFILE: printf(" The system-wide limit on the total number of open files has been reached.\n");
                         break;

            case ENOENT: printf(" The O_CREAT flag was not specified in oflag and no semaphore with this name exists; or, O_CREAT was specified, but name wasn't well formed.  \n");
                         break;

            case ENOMEM: printf(" Insufficient memory.  \n");
                         break;
            default:
                         printf("Unknown error\n");
        }
    }
    sem_unlink(name);
    return ret;
}

void thread_semaphore_destroy(struct Semaphore *sem)
{
    sem_close(sem->semaphore);
}

void thread_semaphore_post(struct Semaphore *sem)
{
    sem_post(sem->semaphore);
}

void thread_semaphore_wait(struct Semaphore *sem)
{
    sem_wait(sem->semaphore);
}

