#pragma once
#include <inttypes.h>

struct TimerSession;
struct TimerSession *timer_session_create(void)
;
void timer_session_set_state(const char *state, struct TimerSession *timer_session)
;
void timer_session_pause(struct TimerSession *timer_session)
;
void timer_session_clear(struct TimerSession *timer_session)
;
void timer_session_summary(void (*timer_session_summary_callback)(const char*, uint64_t, uint64_t, uint64_t),
	struct TimerSession *timer_session)
;
void timer_session_print(struct TimerSession *timer_session)
;
uint64_t timer_session_get_total_ticks(struct TimerSession *timer_session)
;
uint64_t timer_session_get_frequency(struct TimerSession *timer_session)
;

//Global timer session for convenience
void timer_session_set_state_g(const char *state)
;
void timer_session_pause_g(void)
;
void timer_session_clear_g()
;
void timer_session_summary_g(void (*timer_session_summary_callback)(const char*, uint64_t, uint64_t, uint64_t))
;
void timer_session_print_g(void)
;
uint64_t timer_session_get_total_ticks_g(void)
;
//NOTE(Vidar):Returns the number of ticks per second
uint64_t timer_session_get_frequency_g(void)
;
