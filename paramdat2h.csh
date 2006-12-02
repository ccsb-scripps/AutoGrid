#!/bin/csh -f
#
echo 'char *param_string[MAX_LINES] = {'
egrep -v '^#|^$' AD4_parameters.dat | sed 's/\(.*\)$/"\1\\n\\0", /'
switch( `uname -s` )
    case "IRIX64":
    case "CYGWIN_NT-5.1":
        echo '"\\0" };'
        breaksw
    default:
        echo '"\0" };'
        breaksw
endsw
echo '// EOF'
