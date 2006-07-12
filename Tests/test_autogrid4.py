#
# 
#
# $Id: test_autogrid4.py,v 1.5 2006/07/12 17:37:23 rhuey Exp $
#
"""
Test AutoGrid.
"""


import os
import types
import unittest
from string import split, strip

built_maps = False
built_maps_no_parameter_library = False
built_maps_no_receptor_types = False #parse all the receptor types from the receptor input file
built_maps_minus_two_types = False  #parse two extra receptor types from the receptor input file

class Autogrid4_hsg1_sm_test(unittest.TestCase):
    
    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps
        self.autogrid = "../autogrid4"

        if not built_maps:
            # Make sure you remove all the products of AutoGrid from
            # any previous tests.
            command = "rm -f hsg1_sm.*map*"
            os.system(command)
            #print "removed prior maps;",
            gpf_filename = 'hsg1_sm.gpf'
            glg_filename = 'hsg1_sm.glg'
            # run autogrid4
            cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
            #print "compute new maps:\n", cmd_str
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            ###print 'waiting...'
            os.wait() # for the child process to finish
            ###print "after wait"
            built_maps = True


    def tearDown(self):
        #??remove map here
        pass


    def run_cmd(self, cmd_str):
        global built_maps
        (i,o,e) = os.popen3(cmd_str) # trap all the outputs
        #print 'calculating maps...'
        os.wait() # for the child process to finish
        built_maps = True


    def compare_autogrid4_maps(self, stem, maptype, precision=2):
        """Compare autogrid4_maps from C vs python
            allow .225 vs .226 to pass using 'precision'
            allow .526 vs .524 to pass using 0.5 'percentage' check
        """
        # read in 'c' and 'py' maps
        c_autogrid4_map = stem + '.' + maptype +'.map'
        fptr = open(c_autogrid4_map)
        c_autogrid4_lines = fptr.readlines()
        fptr.close()
        #read in python AutoGrid4 results
        py_autogrid4_map = 'py_' + stem +'.'+ maptype+'.map'
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
                self.assertAlmostEquals(c_num, py_num, precision)


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



class Autogrid4_hsg1_sm_no_parameter_library_test(Autogrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps_no_parameter_library
        self.autogrid = "../autogrid4"

        if not built_maps_no_parameter_library:
            # Make sure you remove all the products of AutoGrid from
            # any previous tests.
            command = "rm -f hsg1_sm.*map*"
            os.system(command)
            #print "removed prior maps;",
            gpf_filename = 'hsg1_sm_no_parameter_file.gpf'
            glg_filename = 'hsg1_sm.glg'
            # run autogrid4
            cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
            #print "compute new maps:\n", cmd_str
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            ###print 'waiting...'
            os.wait() # for the child process to finish
            ###print "after wait"
            built_maps_no_parameter_library = True


class Autogrid4_hsg1_sm_no_receptor_types_test(Autogrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps_no_receptor_types
        self.autogrid = "../autogrid4"

        if not built_maps_no_receptor_types:
            # Make sure you remove all the products of AutoGrid from
            # any previous tests.
            command = "rm -f hsg1_sm.*map*"
            os.system(command)
            #print "removed prior maps;",
            gpf_filename = 'hsg1_no_receptor_types.gpf'
            glg_filename = 'hsg1_sm.glg'
            # run autogrid4
            cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
            #print "compute new maps:\n", cmd_str
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            #print 'waiting...'
            os.wait() # for the child process to finish
            #print "after wait"
            built_maps_no_receptor_types = True


class Autogrid4_hsg1_sm_minus_two_types_test(Autogrid4_hsg1_sm_test):

    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        global built_maps_minus_two_types
        self.autogrid = "../autogrid4"

        if not built_maps_minus_two_types:
            # Make sure you remove all the products of AutoGrid from
            # any previous tests.
            command = "rm -f hsg1_sm.*map*"
            os.system(command)
            #print "removed prior maps;",
            gpf_filename = 'hsg1_no_receptor_types.gpf'
            glg_filename = 'hsg1_sm.glg'
            # run autogrid4
            cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
            #print "compute new maps:\n", cmd_str
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            #print 'waiting...'
            os.wait() # for the child process to finish
            #print "after wait"
            built_maps_minus_two_types = True



if __name__ == '__main__':
    test_cases = [
        'Autogrid4_hsg1_sm_test',
        'Autogrid4_hsg1_sm_no_parameter_library_test',
        'Autogrid4_hsg1_sm_no_receptor_types_test',
        'Autogrid4_hsg1_sm_minus_two_types_test',
    ]
    unittest.main( argv=([__name__,] + test_cases))  # non-verbose output
    # optional:  for verbose output, use this:
    # unittest.main( argv=([__name__, '-v'] + test_cases)) # verbose output

