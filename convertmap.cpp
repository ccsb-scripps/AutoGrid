
{
    /****  MANIPULATION CODE  ****/
    char numout;

    if (num == 0.){
	numout = 0;
    } else if ((-12.8 < num) && (num < 0.)) {
	numout = num * 10.;
    } else if ((0. < num) && (num < 1280.)) {
	numout = num / 10.;
    } else if (num >= 1280.)
	numout = 127;
    else
	numout = -128;

    return numout;
}
