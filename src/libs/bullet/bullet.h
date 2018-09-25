#pragma once
#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

EXTERNC void test_bullet(void);

#undef EXTERNC
