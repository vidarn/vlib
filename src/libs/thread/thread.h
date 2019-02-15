#pragma once
/*
#ifdef _WIN32
#else
#include "thread_pthread.h"
#endif
*/

struct Semaphore;
struct ThreadHandle;

struct ThreadHandle *thread_start(unsigned long (*func)(void *),void *param);
void thread_wait(struct ThreadHandle *handle);
struct Semaphore *thread_semaphore_create(unsigned int val, const char *name);
void thread_semaphore_destroy(struct Semaphore *sem);
void thread_semaphore_post(struct Semaphore *sem);
void thread_semaphore_wait(struct Semaphore *sem);


long thread_atomic_increment(long volatile *i);