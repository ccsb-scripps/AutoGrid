#define MAXKEY (256*256)
#include <stdlib.h>
#include <string.h>
#include "structs.h" // needed for ParameterEntry structure
#include "atom_parameter_manager.h"
typedef ParameterEntry PE;
static PE *dictionary[MAXKEY];
static unsigned int hash(const char key[]) {
    switch (strlen(key)) {
        case 0: return 0;
        case 1: return (unsigned int)key[0];
        default: return (unsigned int)key[0] + 256*(unsigned int)key[1];
    }
}
void apm_enter(const char key[], PE value) {
    if (dictionary[hash(key)] == NULL) {
        dictionary[hash(key)] = (PE *) calloc(1, sizeof(PE));
    }
    *(dictionary[hash(key)]) = value;  // this replaces, as well as inserts
    return;
}
PE * apm_find(const char key[]) {
    return dictionary[hash(key)];
}
