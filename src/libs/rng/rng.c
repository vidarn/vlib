#include "pcg_basic.h"
#include "rng.h"
#include <stdlib.h>
#ifdef _MSC_VER
#include <malloc.h>
#endif

struct RNGState{
    pcg32_random_t r;
};


struct RNGState *rng_new_random_state(uint64_t seed)
{
    struct RNGState *ret = calloc(1,sizeof(struct RNGState));
    pcg32_srandom_r(&ret->r, seed, (intptr_t)ret);
    return ret;
}

struct RNGState *rng_clone_state(struct RNGState *state)
{
    struct RNGState *ret = calloc(1,sizeof(struct RNGState));
    ret->r = state->r;
    return ret;
}
    
uint32_t rng_uniform(struct RNGState *state)
{
    return pcg32_random_r(&state->r);
}

uint32_t rng_below(struct RNGState *state, uint32_t max_plus_one)
{
    return pcg32_boundedrand_r(&state->r, max_plus_one);
}

uint32_t rng_range(struct RNGState *state, uint32_t min, uint32_t max_plus_one)
{
    return min + pcg32_boundedrand_r(&state->r, max_plus_one - min);
}

void rng_shuffle(uint8_t *a, uint32_t n, uint32_t element_size,
    struct RNGState *state)
{
    uint8_t *tmp = alloca(element_size);
    for(int i=0;i<n-1;i++){
        int j = rng_range(state, i, n);
        memcpy(tmp,a+i*element_size,element_size);
        memcpy(a+i*element_size,a+j*element_size,element_size);
        memcpy(a+j*element_size,tmp,element_size);
    }
}
