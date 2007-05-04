#
# 
#
# $Id: test_autogrid4.py,v 1.11 2007/05/04 08:00:38 garrett Exp $
#
"""
Test AutoGrid.
"""

import sys
import os
import types
import unittest
import getopt
from string import split, strip

#______________________________________________________________________________
#
# Global variables

autogrid_executable = "../autogrid4" # where the AutoGrid executable resides
gpf_directory = '.' # where the input GPF files reside
test_output_directory = '.' # where the GLG files will be written

try:
    opts, argv = getopt.getopt(sys.argv[1:], "g:e:o:",
    ["gpf-directory=","executable=","test-output-directory="]) 
except getopt.GetoptError, v:
    usage() 
    sys.exit(2)

for o,a in opts:
    if o in ("-d", "--gpf-directory"):
        gpf_directory = a
    if o in ("-e", "--executable"):
        autogrid_executable = a
    if o in ("-o","--test-output-directory"):
        test_output_directory = a


built_maps = False
#no parameter library keyword 
built_maps_no_parameter_library = False
#parse all the receptor types from the receptor input file
built_maps_no_receptor_types = False 
#parse two extra receptor types from the receptor input file
built_maps_minus_two_types = False  
#no receptor types, ligand_types keyword preceeds receptor filename
built_maps_ligand_types_before_receptor = False 

#______________________________________________________________________________

def usage():
    """Print out the usage of this command."""
    print """Usage:  python test_autogrid4.py [-g <string>] [-e <string>] [-o <string>]

where:
    -g, --gpf-directory
        specifies the directory containing the GPFs to be tested;
        this flag is optional; default is '.'
    -e, --executable
        specifies the path to the AutoGrid executable to be tested;
        this flag is optional; default is '../autogrid4'
    -o, --test-output-directory
        specifies the directory where the output GLGs will be written;
        this flag is optional; default is '.'

NOTE:  these may be relative to the directory where this script was invoked.
"""

#______________________________________________________________________________

def run_AutoGrid( gpf_filename, glg_filename ):
    """Launch AutoGrid, using the specified AutoGrid executable and GPF,
    create the specified GLG, and trap all the outputs from standard output
    and standard error."""
    gpf = gpf_directory + os.sep + gpf_filename
    glg = test_output_directory + os.sep + glg_filename
    command = "rm -f " + glg
    os.system( command )
    command = "%s -p %s -l %s" % ( autogrid_executable, gpf, glg )
    print '\nRunning ' + autogrid_executable + ' using GPF "'+gpf+'", saving results in "'+glg+'":'
    try:
        ( i, o, e ) = os.popen3( command ) # trap all the outputs
        os.wait() # for the child process to finish
        return True
    except:
        print "\nUnable to run " + autogrid_executable + "."
        return False

def rm( filename ):
    """Removes files in the test_output_directory."""
    ###os.system( "rm -f "+test_output_directory+"/"+filename )

#______________________________________________________________________________

class AutoGrid4_hsg1_sm_test(unittest.TestCase):
    
    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps
        self.autogrid = autogrid_executable

        if not built_maps:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            rm("hsg1_sm.*map*")
            built_maps = run_AutoGrid('hsg1_sm.gpf', 'hsg1_sm.glg')

    def tearDown(self):
        #??remove map here
        pass

   # def run_cmd(self, cmd_str):
   #   global built_maps
   #     (i,o,e) = os.popen3(cmd_str) # trap all the outputs
   #     #print 'calculating maps...'
   #     os.wait() # for the child process to finish
   #     built_maps = True

    def compare_autogrid4_maps(self, stem, maptype, precision=2):
        """Compare autogrid4_maps from C vs python
            allow .225 vs .226 to pass using 'precision'
            allow .526 vs .524 to pass using 0.5 'percentage' check
        """
        # read in 'c' and 'py' maps
        c_autogrid4_map = test_output_directory + '/' + stem + '.' + maptype +'.map'
        fptr = open(c_autogrid4_map)
        c_autogrid4_lines = fptr.readlines()
        fptr.close()
        #read in python AutoGrid4 results
        py_autogrid4_map = gpf_directory + '/' + 'py_' + stem +'.'+ maptype+'.map'
        fptr = open(py_autogrid4_map)
        py_autogrid4_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(c_autogrid4_lines), len(py_autogrid4_lines))
        for atomL, py_atomL in zip(c_autogrid4_lines[6:], py_autogrid4_lines[6:]):
            c_num = float(strip(atomL))
            py_num = float(strip(py_atomL))
            if abs(c_num)>.30:
                # use percentage difference for large(r) values:
                #  allowing up to 0.5% error 
                #if ((abs(c_num-py_num)/(c_num))*100.<.5)==False:
                #    print  "c_num=", c_num," py_num=", py_num, "test=", (abs(c_num-py_num)/(c_num))*100.
                self.assertEquals(((abs(c_num-py_num)/(c_num))*100.)<.5, True)
            else:
                #but use precision only for smaller values
                cutoff = 10**(-precision)
                self.assertEqual(abs(c_num-py_num)<cutoff, True)
                #self.assertAlmostEquals(c_num, py_num, precision)
        # clean up
        rm(stem + "." + maptype + ".map")

    def test_hsg1_estat(self):
        """test_hsg1 estat case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'e')

    def test_hsg1_dsolv(self):
        """test_hsg1 dsolv case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'd')

    def test_hsg1_A(self):
        """test_hsg1 A case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'A')

    def test_hsg1_C(self):
        """test_hsg1 C case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'C')
 
    def test_hsg1_N(self):
        """test_hsg1 N case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'N')

    def test_hsg1_HD(self):
        """test_hsg1 HD case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'HD')

    def test_hsg1_NA(self):
        """test_hsg1 NA case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'NA')

    def test_hsg1_OA(self):
        """test_hsg1 OA case"""
        # compare resulting map with saved map
        self.compare_autogrid4_maps("hsg1_sm", 'OA')


class AutoGrid4_hsg1_sm_no_parameter_library_test(AutoGrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_no_parameter_library
        self.autogrid = autogrid_executable
        if not built_maps_no_parameter_library:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            rm("hsg1_sm.*map*")
            built_maps_no_parameter_library = run_AutoGrid( 'hsg1_sm_no_parameter_file.gpf', 'hsg1_sm.glg' )


class AutoGrid4_hsg1_sm_no_receptor_types_test(AutoGrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_no_receptor_types
        self.autogrid = autogrid_executable
        if not built_maps_no_receptor_types:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            rm("hsg1_sm.*map*")
            built_maps_no_receptor_types = run_AutoGrid( 'hsg1_no_receptor_types.gpf', 'hsg1_sm.glg' )


class AutoGrid4_hsg1_sm_minus_two_types_test(AutoGrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.  Locate the autogrid binary now during setUp."""
        global built_maps_minus_two_types
        self.autogrid = autogrid_executable
        if not built_maps_minus_two_types:
            # Make sure you remove all the products of AutoGrid from any previous tests.
            rm("hsg1_sm.*map*")
            built_maps_minus_two_types = run_AutoGrid( 'hsg1_no_receptor_types.gpf', 'hsg1_sm.glg' )


#class AutoGrid4_ligand_types_before_receptor_test(AutoGrid4_hsg1_sm_test):

#    def setUp(self):
#        """Set up for autogrid4 tests. Locate the autogrid binary now during setUp."""
#        global built_maps_ligand_types_before_receptor
#        self.autogrid = autogrid_executable
#        if not built_maps_ligand_types_before_receptor:
#            # Make sure you remove all the products of AutoGrid from any previous tests.
#            rm("hsg1_sm.*map*")
#            built_maps_ligand_types_before_receptor = run_AutoGrid( 'hsg1_ligand_types_before_receptor.gpf', 'hsg1_sm.glg')


if __name__ == '__main__':
    test_cases = [
        'AutoGrid4_hsg1_sm_test',
        'AutoGrid4_hsg1_sm_no_parameter_library_test',
        'AutoGrid4_hsg1_sm_no_receptor_types_test',
        'AutoGrid4_hsg1_sm_minus_two_types_test',
        #'AutoGrid4_ligand_types_before_receptor_test',
    ]
    unittest.main( argv=([__name__,] + test_cases))  # non-verbose output
    # optional:  for verbose output, use this:
    # unittest.main( argv=([__name__, '-v'] + test_cases)) # verbose output

# EOF
