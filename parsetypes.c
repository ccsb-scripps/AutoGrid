/* parsetypes.c */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "autogrid.h"

int parsetypes(char * line, char *words[], int maxwords)
/*utility func for parsing types*/
{
/******************************************************************************/
/*      Name: parsetypes                                                      */
/*  Function: Parse the AutoGrid types lines                                  */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 02/01/95 (1-feb-1995)                                           */
/*----------------------------------------------------------------------------*/
/*    Inputs: line, array of pointers, cut-off number of words                */
/*   Returns: integer, number of types found.                                 */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/02/03 RH      Entered code.                                             */
/******************************************************************************/

    char *p = line;
    int nwords = 0;
    int found_keyword = 0;
    int index = 0;

    while(1) {
        while(isspace(*p)){
            p++;
            index++;
        };
        if (*p == '\0'){
            return nwords;
        };
        /*words[nwords++] = p;*/
        if(found_keyword==0){
            found_keyword++;
        } else {
            words[nwords++] = p;
        };
        while(!isspace(*p) && *p!='\0'){
            p++;
            index++;
        };
        if(*p=='\0'){
            return nwords;
        };
        *p++ = '\0';
        index++;
        if(nwords >=maxwords){
            return nwords;
        };
    }

}
/* EOF */
