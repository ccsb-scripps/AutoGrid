#!/bin/csh -f
#
echo 'char *param_string[MAX_LINES] = {'
egrep -v '^#|^$' AD4_parameters.dat | sed 's/\(.*\)$/"\1\\n\\0", /'
echo '"\0" };'
echo '// EOF'
