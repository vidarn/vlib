#pragma once

struct RNGState;
#include <inttypes.h>

struct RNGState *rng_new_random_state(uint64_t seed);
struct RNGState *rng_clone_state(struct RNGState *state);
uint32_t rng_uniform(struct RNGState *state);
uint32_t rng_below(struct RNGState *state, uint32_t max_plus_one);
uint32_t rng_range(struct RNGState *state, uint32_t min, uint32_t max_plus_one);

void rng_shuffle(uint8_t *a, uint32_t n, uint32_t element_size,
    struct RNGState *state);
