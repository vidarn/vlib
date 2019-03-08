#include "timer.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef WIN32
#include <Windows.h>
static uint64_t timer_get_tick()
{
	LARGE_INTEGER l;
	QueryPerformanceCounter(&l);
	return l.QuadPart;
}
// The current performance-counter frequency, in counts per second
static uint64_t timer_get_frequency()
{
	LARGE_INTEGER l;
	QueryPerformanceFrequency(&l);
	return l.QuadPart;
}
#else
#include <time.h>
static uint64_t timer_get_tick()
{
	struct timespec time;
	clock_gettime(CLOCK_REALTIME, &time);
	return time.tv_sec * (uint64_t)1000000000L + time.tv_nsec;
}
static uint64_t timer_get_frequency()
{
	return (uint64_t)1000000000L ;
}
#endif


struct TimerSession {
	const char **state_names;
	uint64_t *state_ticks;
	int alloc_states;
	int num_states;
	uint64_t last_tick;
	int current_state;
};

struct TimerSession *timer_session_create(void)
{
	struct TimerSession *ret = calloc(1, sizeof(struct TimerSession));
	return ret;
}

void timer_session_set_state(const char *state, struct TimerSession *timer_session)
{
	uint64_t tick = timer_get_tick();

	if (timer_session->alloc_states == 0) {
		timer_session->alloc_states = 128;
		timer_session->state_names = calloc(timer_session->alloc_states, sizeof(const char *));
		timer_session->state_ticks = calloc(timer_session->alloc_states, sizeof(int));
		timer_session->current_state = -1;
	}

	int last_state_index = timer_session->current_state;
	if (last_state_index >= 0) {
		timer_session->state_ticks[last_state_index] += tick - timer_session->last_tick;
	}

	int state_index = -1;
	if (state) {
		for (int i = 0; i < timer_session->num_states; i++) {
			if (strcmp(timer_session->state_names[i], state) == 0) {
				state_index = i;
			}
		}
		if (state_index == -1) {
			state_index = timer_session->num_states++;
			if (timer_session->num_states > timer_session->alloc_states) {
				timer_session->alloc_states *= 2;
				timer_session->state_names = realloc(timer_session->state_names, timer_session->alloc_states* sizeof(const char*));
				timer_session->state_ticks = realloc(timer_session->state_ticks, timer_session->alloc_states* sizeof(int));
			}
			timer_session->state_names[state_index] = strdup(state);
			timer_session->state_ticks[state_index] = 0;
		}
	}
	timer_session->current_state = state_index;
	
	timer_session->last_tick = timer_get_tick();
}

void timer_session_pause(struct TimerSession *timer_session)
{
	timer_session_set_state(0, timer_session);
}

void timer_session_summary(void (*timer_session_summary_callback)(const char*, uint64_t, uint64_t, uint64_t),
	struct TimerSession *timer_session)
{
	timer_session_pause(timer_session);
	uint64_t total_ticks = 0;
	for (int i = 0; i < timer_session->num_states; i++) {
		total_ticks += timer_session->state_ticks[i];
	}
	uint64_t freq = timer_get_frequency();
	for (int i = 0; i < timer_session->num_states; i++) {
		timer_session_summary_callback(timer_session->state_names[i], timer_session->state_ticks[i],
			freq, total_ticks);
	}
}

static void timer_session_summary_print_callback(const char *state, uint64_t ticks, uint64_t frequency, uint64_t total_ticks)
{
	float percent = (float)(100.0*(double)ticks / (double)total_ticks);
	printf("%s: %g%%\n", state, percent);
}

void timer_session_print(struct TimerSession *timer_session)
{
	timer_session_summary(timer_session_summary_print_callback, timer_session);
}

static struct TimerSession g_timer_session = { 0 };

void timer_session_set_state_g(const char *state)
{
	timer_session_set_state(state, &g_timer_session);
}
void timer_session_pause_g()
{
	timer_session_pause(&g_timer_session);
}
void timer_session_summary_g(void(*timer_session_summary_callback)(const char*, uint64_t, uint64_t, uint64_t))
{
	timer_session_summary(timer_session_summary_callback, &g_timer_session);
}
void timer_session_print_g()
{
	timer_session_print(&g_timer_session);
}
