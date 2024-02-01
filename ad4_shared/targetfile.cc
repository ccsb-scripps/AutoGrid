/* 
 $Id: targetfile.cc,v 1.1 2020/05/21 15:26:38 mp Exp $
 
 Support for "target" files as input to AutoDock; a target file is a zipped
 collection of map files and other information.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*******************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

static const char* const ident[] = {ident[1], "@(#)$Id: targetfile.cc,v 1.1 2020/05/21 15:26:38 mp Exp $"};

#include <stdlib.h>
#include <stdio.h> // for fgets()
#include <ctype.h> // for isspace()
#include <string.h>
#include <unistd.h> // for mkstemps()
#include "autocomm.h"
#include "parse_dpf_line.h"
#include "dpftoken.h"
#include "constants.h" // for LOGVV.. defs
#include "stop.h"
#include "targetfile.h"

#ifndef HAVE_ZIP_H
	// provide functions to indicate target files are not supported
	//
int target_file_open ( FILE * dpf_in, FILE ** dpf_out, char *dpf_out_filename, 
	const int outlev, FILE *logFile ) {
	return 0; // error
}
int target_file_remove_tfiles ( const char *arg) {
	return 0; // OK
}
int target_file_capability() {
	return 0; // no capability
}
#else

#include <zip.h> // for the libzip functions

// #define DEBUGMAIN

#define streq(a,b) (0==strcasecmp(a,b)) // case-independent string match

// convenience macro for parsing the (many) single-argument DPF lines:
//  if not 1 (non-ignored) argument, stop, reporting fatal error
#define get1arg(line, fmt, addr, token) if(1!=sscanf(line, fmt, addr))stop("syntax error in " token " line")
// convenience macro for reporting syntax errors in DPF lines:
#define syntaxstop(s) {char ss[LINE_LEN+50];sprintf(ss,"syntax error or illegal value in %s line",s);stop(ss);}
// convenience macro for looking at ends of file names (case-sensitive)
#define str_ends_in(s, ending) (0==strcmp(ending, s+strlen(s)-strlen(ending)))

// convenience macro for plural noun string
#define pl(i) ((i==1)?"":"s")

#include <sys/param.h>

/* globals : */
extern int debug;


#ifdef DEBUGMAIN
int debug;
#define stop(s) fprintf(stderr, "%s", s),exit(-1)
int
main(int argc,char **argv) {
char * dpf = "t2.dpf";
if (argc==2) {
	dpf=argv[2];
}

FILE *dpfp;
FILE *newdpfp = NULL;
char newdpfname[PATH_MAX];
dpfp = fopen(dpf,"r");
if (NULL == dpfp) {
	fprintf(stderr,"Cant open old dpf %s\n", dpf);
	exit (1);
	}
int rc;
int outlev=4;

rc = target_file_open( dpfp, &newdpfp, newdpfname, outlev, stderr);
char line[1000];

if(debug)fprintf(stderr," target_file_open rc=%d newdpfname=%s newdfpp=%p\n", 
  rc, newdpfname, newdpfp);

// read back the new DPF - note no rewind needed here
while ( fgets(line, 1000, newdpfp)!=NULL ) {
	fprintf(stderr, "%s", line);
	}
printf(" calling target_file_remove_tfiles \n");
rc = target_file_remove_tfiles("all");
printf(" target_file_remove_tfiles rc=%d\n", rc);
}
#endif

//
//  zip archive file entries are extracted and copied to file system temp space
//  for AutoDock to use. This includes .map, .fld, and potentially .dat types.
//  A new (substitute) Docking Parameter File (PDF) is created and a pointer
//   to the open file is returned to AutoDock main().
//  AutoDock main should call  target_file_remove_tfiles("all") after all
//   the files are read. The code below inserts this command into the new DPF.
//
// The "character" type maps [a byte-binary format] are not supported 
//   because I do not think they are used or even working in AutoDock.
//
//  M Pique, May 2020


// filename max length is taken from system include file
// MAX_MAPS is an AutoDock compile-time constant (autocomm.h):
//    MAX_ATOM_TYPES + NUM_NON_VDW_MAPS + MAX_MAPS_PAD 
// NUM_NON_VDW_MAPS is 2: electrostatic and desolvation
// PAD is currently zero)
// We add 2 here for the fld file and a future .dat file
#define MAX_NTFILES (MAX_MAPS+2)




// variables shared within target_file.cc 
static char *tfile_name[MAX_NTFILES]; // map name as given in the target file
static char *tfile_path[MAX_NTFILES]; // path to each uncompressed tfile
static char map_type_name[MAX_NTFILES][3]; // atom type names for each map (not fld)
static int ntfiles; // includes all .map, incl  e,d maps if any, and fld

static char dpf_out_name[PATH_MAX];
static FILE * dpf_out_fptr;

// declarations for local helper functions, see end of this source file:
static int target_file_generate ( const char *target_file_name, 
  FILE * dpf_fptr,  // must be open for writing
  const int outlev, FILE *logFile);
static int get_map_type_str (char map_name[], char *map_type_name /* char[3] */ ) ;

int
target_file_open ( FILE * dpf_in, FILE ** dpf_out, char *dpf_out_filename, 
	const int outlev, FILE *logFile ) 
{
/* copy dpf_in to (new) dpf_out, modifying TARGET line, set ptr and name of new file,
 * return 1 if OK, 0 if fatal error  */

char line[LINE_LEN];
char target_file_name[LINE_LEN];

int dpf_keyword = -1;
int generate_dpf_rc = -1;


/* create new edited dpf */
(void) strcpy(dpf_out_name, "/tmp/AutoDock.XXXXXX.dpf");
mkstemps ( dpf_out_name, strlen(".dpf"));
dpf_out_fptr = fopen (dpf_out_name, "w+");
if (NULL==dpf_out_fptr) {
	fprintf(logFile, "target_file open could not create new dpf %s\n", dpf_out_name);
	return 0; // fatal
}
	

/* read incoming DPF - most lines are passed unchanged, some trigger warning,
   "target" triggers action.
*/

while( fgets(line, LINE_LEN, dpf_in) != NULL ) { 
#ifdef DEBUGMAIN
//fprintf(stderr, "not parsing dpf line %s\n",line);fflush(stderr); 
int nmatch;
char token[30], value[30];
nmatch=sscanf(line, "%s %s", token, target_file_name);
if (nmatch==2&&streq("target",token)) dpf_keyword=DPF_TARGETFILE;
else dpf_keyword=0; // other

#else
    dpf_keyword = parse_dpf_line( line );
#endif
    if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0';  // remove newline if last char in line

    switch( dpf_keyword ) {
    case DPF_FLD:
    case DPF_MAP:
    case DPF_ELECMAP:
    case DPF_DESOLVMAP:
    case DPF_LIGAND_TYPES:
#define DATFILE
#ifdef DATFILE
    case DPF_PARAMETER_LIBRARY:
#endif
	
	if (outlev>LOGMIN) {
		fprintf (logFile, "Warning, ignoring '%s' as superceded by 'target'\n",
			line);
		}
	break;

    case DPF_TARGETFILE:

        get1arg( line, "%*s %s", target_file_name, "TARGET" );
	generate_dpf_rc = target_file_generate ( target_file_name, dpf_out_fptr, outlev, logFile);
        break;
    default:
        fprintf (dpf_out_fptr, "%s\n", line);
	break;
    } // switch( dpf_keyword )
} // while

// Rewind new DPF, so caller can start reading it
rewind (dpf_out_fptr);


// if fatal error returned by target_file_generate(), delete all
// the temporary files it made and remove the new DPF.
if(generate_dpf_rc==0) {
	unlink(dpf_out_name);
	target_file_remove_tfiles("all");
	return 0; // fatal
}


// set return values:
*dpf_out = dpf_out_fptr;
strncpy(dpf_out_filename, dpf_out_name, PATH_MAX);

// remove our new DPF - it will stay on disk until caller does fclose or exit
unlink (dpf_out_name);
return 1; // OK
} // target_file_open function


int
target_file_generate ( const char *target_file_name, 
  FILE * dpf_fptr,  // must be open for writing
  const int outlev, FILE *logFile)
{
static Boole B_have_maps = FALSE;
static Boole B_have_fld = FALSE;
static Boole B_have_dat = FALSE;


struct zip *zip;
int ziperr;


zip = zip_open (target_file_name, 0, &ziperr);

if (zip==NULL) {
        fprintf(logFile, "target file \"%s\" zipopen fail %d\n", target_file_name, ziperr);
        return 0; // fatal
        }
//int64_t nentries;

zip_int64_t nentries;
//nentries = zip_get_num_entries (zip, (zip_flags_t) ZIP_FL_UNCHANGED);
nentries = zip_get_num_entries (zip, 0);
//fprintf(stderr, "target %s %lld entries\n", target_file_name, nentries);

int i;
/* debug
for(i=0;i<nentries;i++) {
	fprintf(stderr, "%2d %s entries\n", 
	i, zip_get_name (zip, (zip_uint64_t) i, 0));
	//i, zip_get_name (zip, (zip_uint64_t) i, (zip_flags_t)0));
	}
*/


for(i=0;i<nentries;i++) {
	//char *ename =  zip_get_name (zip, (zip_uint64_t) i, (zip_flags_t)0);
	char *ename =  (char *) zip_get_name (zip, (zip_uint64_t) i, 0);


	if (! (str_ends_in(ename, ".fld") 
#ifdef DATFILE
		|| str_ends_in(ename, ".dat")
#endif
		|| str_ends_in(ename, ".map"))) continue;
		// extract map and fld files from archive, copy to new location

		if (ntfiles>=MAX_NTFILES) {
			fprintf(logFile, "target file %s has too many maps", target_file_name);
			return 0;  // fatal
			}

		//zip_file_t *zf = zip_fopen_index (zip, (zip_uint64_t) i, (zip_flags_t) 0);
		struct zip_file *zf = zip_fopen_index (zip, (zip_uint64_t) i,  0);
		if (str_ends_in(ename, ".map")) {
			int rc;
			rc = get_map_type_str(ename, map_type_name[ntfiles]);
			if(rc==0) { // error
				fprintf(logFile, "Illegal map type %s in entry %s from target file %s",
					map_type_name[ntfiles], ename, target_file_name);
				return 0; // fatal
				}
			B_have_maps = TRUE;
			}

		// now copy the entry to a temporary file for AutoDock to read
		char suffix[20];
		if (str_ends_in(ename,".map")) sprintf(suffix, ".%s.map", map_type_name[ntfiles]);
		else if (str_ends_in(ename,".fld")) sprintf(suffix, ".fld");
#ifdef DATFILE
		else if (str_ends_in(ename,".dat")) sprintf(suffix, ".dat");
#endif
		else sprintf(suffix, "%s", "" );  // can't happen. Guard is above...

		char tempfilename[PATH_MAX];
		int tempfiledesc = -1;
		//sprintf (tempfilename, "%s/AutoDockTargetXXXXXXX%s",
			//getenv("TMPDIR")!=NULL?getenv("TMPDIR"):"/tmp", suffix);
		sprintf (tempfilename, "/tmp/AutoDockTargetXXXXXXX%s",
			suffix);
		if(outlev>LOGRUNVVV) fprintf(stderr, "tempfilename pre  = '%s' fd=%d\n", 
			tempfilename, tempfiledesc);
		tempfiledesc = mkstemps( tempfilename, strlen(suffix));
		if (tempfiledesc<0) {
			fprintf( logFile, "target file could not create temporary file \"%s\" rc=%d\n",
				tempfilename, tempfiledesc);
			return 0; // fatal
		}
		if(outlev>LOGRUNVVV) fprintf(stderr, "tempfilename post = '%s' fd=%d\n", 
		  tempfilename, tempfiledesc);

		// save these names into our local static storage
		tfile_name[ntfiles] = (char *) malloc( 1 + strlen(ename)); // TODO check malloc
		strncpy(tfile_name[ntfiles],ename,1+strlen(ename));
		tfile_path[ntfiles] = (char *) malloc( 1 + strlen(tempfilename)); // TODO check malloc
		strncpy(tfile_path[ntfiles],tempfilename,1+strlen(tempfilename));


#define BUFSIZE (8*512*512)
		char buf[BUFSIZE]; // for copying to temp file
		int nread;
		while ((nread = zip_fread (zf, buf, (zip_uint64_t) BUFSIZE))>0) {
			int nwrite = write (tempfiledesc, buf, nread);
			if (nwrite != nread ) {
				fprintf(logFile ,"target expansion error copying to tmp file");
				return 0; // fatal
			}
			if (outlev>LOGRUNVVV) fprintf(logFile, "zip_fread %d bytes, wrote tempfile %d\n", nread, nwrite);
			}
		(void) close (tempfiledesc); // TODO check rc
		(void) zip_fclose (zf); // TODO check rc

/*
		// test temp file:
		FILE * tempfileptr = fopen(tempfilename, "r");
		bzero(buf, BUFSIZE);
		fread (buf, 1, 200, tempfileptr);
		fprintf (stderr, "%s\n", buf);
		fclose(tempfileptr);
*/

#ifdef DATFILE
		if ( str_ends_in(ename, ".dat")) {
			// TODO warn if multiple
			B_have_dat = TRUE;
			}
#endif
		if ( str_ends_in(ename, ".fld")) {
			// TODO warn if multiple
			B_have_fld = TRUE;
			}

		ntfiles++;
		} // next zip entry
	// close the zip file itself, without deleting it or saving changes (none!)
	zip_close(zip);

	// check that we have a fld file and at least one map file
	// if not, report error so caller can remove all the temp files and the DPF under construction
	if (! (B_have_fld && B_have_maps
#ifdef DATFILE
	&& B_have_dat
#endif
)) {
		if (!B_have_fld) fprintf( logFile,  "no fld file in target file \"%s\"\n", target_file_name);
		if (!B_have_maps) fprintf( logFile, "no map files in target file \"%s\"\n", target_file_name);
#ifdef DATFILE
		if (!B_have_dat) fprintf( logFile, "no dat file in target file \"%s\"\n", target_file_name);
#endif
		return 0; // fatal
	}
#ifdef DATFILE
	// generate dat line in new dpf
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".dat")) {
		fprintf (dpf_fptr, 
		 "parameter_file %s\t# fld in %s from target file %s\n", 
		   tfile_path[m], tfile_name[m], target_file_name);
		}
#endif

	 /*
	 * generate ligand types line in new dpf (do not include d or e map)
         *  ligand_types C HD OA P               # ligand atom type names
         *
         *  The order of the arguments is the index that AutoDock
         *  will use for look up in the grid maps, "map_index".
         */
	fprintf (dpf_fptr, "ligand_types ");
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".map")) {
		if (streq(map_type_name[m],"d") ||
		 streq(map_type_name[m],"e")) continue;
		fprintf (dpf_fptr, " %s", map_type_name[m]);
		}
	fprintf (dpf_fptr, "\t\t# from target file %s\n", target_file_name);

	// generate fld line in new dpf
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".fld")) {
		fprintf (dpf_fptr, 
		 "fld %s\t# fld in %s from target file %s\n", 
		   tfile_path[m], tfile_name[m], target_file_name);
		}

	// generate map lines in new dpf - first all 'non special'atom type lines, 
	//   then 'e' then 'd' maps. TODO check for duplicate c,d,or e
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".map")) {
		if (streq(map_type_name[m],"d") ||
		 streq(map_type_name[m],"e")) continue;
		fprintf (dpf_fptr, 
		 "map %s\t# %2s in %s from target file %s\n", 
		   tfile_path[m], map_type_name[m], tfile_name[m], target_file_name);
		}
	// elecmap if present
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".map") &&
	    streq(map_type_name[m],"e")) {
		fprintf (dpf_fptr, 
		 "elecmap %s\t# %2s in %s from target file %s\n", 
		   tfile_path[m], map_type_name[m], tfile_name[m], target_file_name);
		}
	// desolvmap if present
	for (int m=0;m<ntfiles;m++) if (str_ends_in(tfile_name[m],".map") &&
	    streq(map_type_name[m],"d")) {
		fprintf (dpf_fptr, 
		 "desolvmap %s\t# %2s in %s from target file %s\n", 
		   tfile_path[m], map_type_name[m], tfile_name[m], target_file_name);
		}
	// Tell AutoDock main to delete the temporary files
	fprintf (dpf_fptr, "target_close all\n");


	
fflush(dpf_fptr);
return 1; // OK
}


int
target_file_remove_tfiles ( const char *arg)
{
/* 
 * Delete a path-specified temporary file or "all"
 * Note this function disregards 'ntfiles' to avoid the complexity of keeping
 *  the tfile_path[] compacted.
 *
 *  It is harmless to call even if no files were opened 
 */
int rc=0; // OK
	for(int m=0;m<MAX_NTFILES;m++) {
		if (tfile_path[m] == NULL) continue;
		if ( ! ( streq(arg,"all") || streq(arg,tfile_path[m]))) continue;
		int r = unlink (tfile_path[m]);
#ifdef NO
		if (outlev>LOGRUNVV || r!=0) {
			fprintf(logFile, "target_file_remove_tfiles unlinking tfile %s, rc=%d\n",
			  map_path[m], r);
			}
#endif
		// free storage and clear static variables:
		free (tfile_name[m]);
		tfile_name[m] = NULL;
		free (tfile_path[m]);
		tfile_path[m] = NULL;
		ntfiles--;
		rc += r;
		}
	return rc;
	}

int
target_file_capability()
{
/*
 * return 0 if this compilation/build cannot read target files,
 * return >0 if it can
 */
return 1;
}
static int
get_map_type_str (char map_name[], char *map_type_name /* char[3] */ ) {
// Extract <begin>*\(<char><optionalchar>\).map into map_type_name, 
//  where <begin> is beginning of string, '.', or '/', and
//  where <char> is not '.', '/', or ' '.
// Minimal legal map names would be  a.map  (of type a)
// or ab.map  (of type ab)
// Not worth setting up a regexp but it would be like
//    [/\.^]*[^/\.[:space:]][^/\.[:space:]]?\.map$
//      START^                            ^END
//
//  return 1 if OK, 0 if error.  If error, put empty string into map_type_name.
//  M Pique, May 2020


map_type_name[0] = '\0';
if (! str_ends_in(map_name, ".map")) return 0; // error
if ( strlen(map_name) < 5) return 0; // error (min would be "a.map")

int typelen = 0;
int bgntype; /* becomes 0-origin index of beginning of the type characters */

for (bgntype=strlen(map_name)-4;  // ie, the character before the final '.' 
     bgntype > 0 
	&& map_name[bgntype-1] != '.'
	&& map_name[bgntype-1] != '/'
	&& (!isspace(map_name[bgntype-1]))
	; typelen++, bgntype--) ;



if(typelen<1||typelen>2) {
	//fprintf(stderr, "bgntype=%d   typelen=%d not 1 or 2\n", bgntype, typelen);
	return 0; // error
	}
if(0) fprintf(stderr, "bgntype=%d '%c'  typelen=%d   '%c'\n",
	bgntype, map_name[bgntype], 
	typelen, map_name[bgntype+typelen]);
int p;
for (p=0;p<typelen;p++) map_type_name[p] = map_name[bgntype+p];
map_type_name[typelen] = '\0';
return 1; // OK

} // end get_map_type_str()

// endif for HAVE_ZIP_H
#endif

#ifdef TESTMAPTYPESTR
main (int argc, char **argv) {
int get_map_type_str (char map_name[], char *map_type_name /* char[3] */ );
/* some maptype test cases: */
int rc;
char *s, t[3];
// OK ones:
printf("should be OK:\n");
s="a.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="Aa.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="s.Ca.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s=".Ca.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a.c..s.ba.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a/Aa.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a/.Aa.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a/a.Aa.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");

// illegal ones:
printf("should be ERROR:\n");
s="";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s=".";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s=".map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="..map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="Aa..map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="s.Abc.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="/.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s=". .map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a. .map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a .map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s=" .map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="./.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a./.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
s="a/.map";rc=get_map_type_str (s, t); printf("%s %s rc=%s\n\n", s,t,rc==0?"ERROR":"OK");
}
#endif
