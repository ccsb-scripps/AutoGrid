
#include <time.h>

struct tms
{
	clock_t tms_utime;	// CPU time used in executing the instructions of the calling process. 
	clock_t tms_stime;	// CPU time used by the system on behalf of the calling process. 
	clock_t tms_cutime;	// sum of the tms_utime values and the tms_cutime values of all terminated child processes of the calling process, whose status has been reported to the parent process by wait or waitpid; see section Process Completion. In other words, it represents the total CPU time used in executing the instructions of all the terminated child processes of the calling process, excluding child processes which have not yet been reported by wait or waitpid. 
	clock_t tms_cstime;	// similar to tms_cutime, but represents the total CPU time used by the system on behalf of all the terminated child processes of the calling process. 
	// All of the times are given in clock ticks. These are absolute values; in a newly created process, they are all zero. See section Creating a Process. 
};


extern clock_t times (struct tms *buffer);

// end