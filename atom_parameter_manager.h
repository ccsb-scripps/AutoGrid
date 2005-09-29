
#ifndef _APM_PROTOTYPES
#   define _APM_PROTOTYPES

#   include "structs.h"
void apm_enter(const char key[], ParameterEntry value);
ParameterEntry * apm_find(const char key[]);

#endif
