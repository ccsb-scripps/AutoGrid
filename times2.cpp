
// times2.c
//
// emulate UNIX times() function for CPU time

#ifdef _WIN32
#include <time.h>

clock_t times (struct tms *buffer)
{
	return clock();
}
//The times function stores the processor time information for the calling process in buffer.
//The return value is the same as the value of clock(): the elapsed real time relative to an arbitrary base. The base is a constant within a particular process, and typically represents the time since system start-up. A value of (clock_t)(-1) is returned to indicate failure.
//Portability Note: The clock function described in section Basic CPU Time Inquiry

#endif   // _WIN32
