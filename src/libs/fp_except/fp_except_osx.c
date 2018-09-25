
#include <fenv.h>
#define ACTIVE_FPE FE_DIVBYZERO| FE_INVALID| FE_OVERFLOW //| FE_UNDERFLOW
#define ENABLE_FPE  feenableexcept (ACTIVE_FPE)
#define DISABLE_FPE fedisableexcept(ACTIVE_FPE)

// Public domain polyfill for feenableexcept on OS X
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c

static inline int feenableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);

    return fesetenv(&fenv) ? -1 : old_excepts;
}

static inline int fedisableexcept(unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    // all previous masks
    unsigned int old_excepts;

    if (fegetenv(&fenv)) {
        return -1;
    }
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= new_excepts << 7;

    return fesetenv(&fenv) ? -1 : old_excepts;
}

void fp_except_enable(void)
{
    ENABLE_FPE;
}

void fp_except_disable(void)
{
    DISABLE_FPE;
}
