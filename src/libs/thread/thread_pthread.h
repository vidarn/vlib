#pragma once
#include <pthread.h>
#include <semaphore.h>
struct ThreadHandle{
    pthread_t handle;
};
struct Semaphore{
    sem_t *semaphore;
};
