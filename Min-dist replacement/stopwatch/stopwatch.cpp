// stopwatch.c
// recording when function starts and stops

#include <stdio.h>
#include <sys/timeb.h>

#include "stopwatch.h"

static struct _timeb start_time, finish_time;

// start the timer
void start(void) {
	_ftime_s(&start_time);
}

// stop the timer and print the time elasped since last start() call
void stop(void) {
	_ftime_s(&finish_time);

	double duration = finish_time.time + 0.001 * finish_time.millitm - (start_time.time + 0.001 * start_time.millitm);
	printf("Time cost = %4.2f\n",duration);
}
