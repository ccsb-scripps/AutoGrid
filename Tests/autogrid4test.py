#
# 
#
# $Id: autogrid4test.py,v 1.5 2004/07/13 17:37:41 rhuey Exp $
#
"""

"""


import os
import types
import unittest
from string import split, strip


class AutoGridTestError(Exception):
    pass


class Autogrid4TestCase(unittest.TestCase):
    
    err_epsilon = 0.006
    #err_epsilon = 0.0005
    #err_epsilon = 0.167  # to compare .005 with .006 FIX THIS!!!
    def assertFloatEquals(self, float1, float2, digits=None):
        if digits is not None:
            float1 = round(float1, digits)
            float2 = round(float2, digits)
        #check whether numbers diff by .001
        #for example .005 and .006
        delta = abs(abs(float1)-abs(float2))
        if abs(float2-float1)==0.0:
            self.assertEquals(1,1)
        elif abs(abs(delta/.001)-1.0) <self.err_epsilon:
            self.assertEquals(1,1)
        else:
            try:
                eps = abs((float1/float2) - 1.0)
            except ZeroDivisionError:
                eps = abs((float2/float1) - 1.0)
            self.assertEquals((eps < self.err_epsilon), True,
                              msg="%f != %f; eps=%f" % (float1, float2, eps))
    # 1 method, 2 names
    assertFloatEqual = assertFloatEquals

    
    def setUp(self):
        """Set up for autogrid4 tests.
        Locate the autogrid binary now during setUp.
        """
        self.autogrid = "../autogrid4"
        # Make sure you remove all the products of AutoGrid from
        # any previous tests.
        command = "rm -f *map*"
        os.system(command)


    def tearDown(self):
        #??remove map here
        pass


    def run_cmd(self, cmd_str):
        (i,o,e) = os.popen3(cmd_str) # trap all the outputs
        #print 'calculating map...'
        os.wait() # for the child process to finish


    def compare_maps(self, stem):
        """Compare maps 
        """
        # read in the old and new maps
        atom_map = 'test' + stem + '.' + stem +'.map'
        estat_map = 'test' + stem + '.e.map'
        dsolv_map = 'test' + stem + '.d.map'
        fptr = open(atom_map)
        atom_lines = fptr.readlines()
        fptr.close()
        fptr = open(estat_map)
        estat_lines = fptr.readlines()
        fptr.close()
        fptr = open(dsolv_map)
        dsolv_lines = fptr.readlines()
        fptr.close()
        #read in the results of the evaluator
        evaluator_filename = 'test' + stem + "_summary.txt"
        fptr = open(evaluator_filename)
        evaluator_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(atom_lines), len(evaluator_lines)+5)
        for atomL, estatL, dsolvL, evalL  in zip(atom_lines[6:], estat_lines[6:],
                                                 dsolv_lines[6:],
                                                 evaluator_lines[1:]):
            #print atomL, 
            eval_list = split(evalL)
            #print eval_list
            self.assertFloatEquals(float(strip(estatL)), round(float(strip(eval_list[0])), 3))
            self.assertFloatEquals(float(strip(atomL)), round(float(strip(eval_list[1])), 3))
            self.assertFloatEquals(float(strip(dsolvL)), round(float(strip(eval_list[2])), 3))


    def compare_e_maps(self, stem):
        """Compare estat maps 
        """
        # read in the old and new maps
        #atom_map = 'test' + stem + '.' + stem +'.map'
        estat_map = 'test' + stem + '.e.map'
        #dsolv_map = 'test' + stem + '.d.map'
        #fptr = open(atom_map)
        #atom_lines = fptr.readlines()
        #fptr.close()
        fptr = open(estat_map)
        estat_lines = fptr.readlines()
        fptr.close()
        #fptr = open(dsolv_map)
        #dsolv_lines = fptr.readlines()
        #fptr.close()
        #read in the results of the evaluator
        evaluator_filename = 'test' + stem + "_summary.txt"
        fptr = open(evaluator_filename)
        evaluator_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(estat_lines), len(evaluator_lines)+5)
        for  estatL, evalL  in zip(estat_lines[6:], evaluator_lines[1:]):
            #print atomL, 
            eval_list = split(evalL)
            #print eval_list
            self.assertFloatEquals(float(strip(estatL)), round(float(strip(eval_list[0])), 3))


    def compare_vdw_maps(self, stem):
        """Compare vdw maps 
        """
        # read in the old and new maps
        atom_map = 'test' + stem + '.' + stem +'.map'
        #estat_map = 'test' + stem + '.e.map'
        #dsolv_map = 'test' + stem + '.d.map'
        fptr = open(atom_map)
        atom_lines = fptr.readlines()
        fptr.close()
        #fptr = open(estat_map)
        #estat_lines = fptr.readlines()
        #ptr.close()
        #fptr = open(dsolv_map)
        #dsolv_lines = fptr.readlines()
        #fptr.close()
        #read in the results of the evaluator
        evaluator_filename = 'test' + stem + "_summary.txt"
        fptr = open(evaluator_filename)
        evaluator_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(atom_lines), len(evaluator_lines)+5)
        ctr = 0
        for  i, x in enumerate(zip(atom_lines[6:], evaluator_lines[1:])):
            atomL = x[0]
            evalL = x[1]
            #print dsolvL, 
            eval_list = split(evalL)
            #print eval_list
            float1 = round(float(strip(atomL)),3)
            float2 = round(float(strip(eval_list[1])),3)
            delta = abs(abs(float1)-abs(float2))
            ok = False
            if abs(float2-float1)==0.0:
                ok = True
            elif abs(abs(delta/.001)-1.0) <self.err_epsilon:
                ok = True
            else:
                try:
                    ok = abs((float1/float2) - 1.0)<self.err_epsilon
                except ZeroDivisionError:
                    ok = abs((float2/float1) - 1.0)<self.err_epsilon
                
            # 1 method, 2 names
            if not ok:     
                print 'i:',i, 'float1=', float1, ':float2=', float2
                ctr = ctr+1
            #self.assertFloatEquals(float(strip(atomL)), round(float(strip(eval_list[1])), 3))
        if ctr!=0: print "ctr=", ctr
        self.assertEquals(ctr, 0)


    def compare_d_maps(self, stem):
        """Compare dsolv maps 
        """
        # read in the old and new maps
        #atom_map = 'test' + stem + '.' + stem +'.map'
        #estat_map = 'test' + stem + '.e.map'
        dsolv_map = 'test' + stem + '.d.map'
        #fptr = open(atom_map)
        #atom_lines = fptr.readlines()
        #fptr.close()
        #fptr = open(estat_map)
        #estat_lines = fptr.readlines()
        #ptr.close()
        fptr = open(dsolv_map)
        dsolv_lines = fptr.readlines()
        fptr.close()
        #read in the results of the evaluator
        evaluator_filename = 'test' + stem + "_summary.txt"
        fptr = open(evaluator_filename)
        evaluator_lines = fptr.readlines()
        fptr.close()
        # skip starting lines
        self.assertEqual(len(dsolv_lines), len(evaluator_lines)+5)
        for  dsolvL, evalL  in zip(dsolv_lines[6:], evaluator_lines[1:]):
            #print dsolvL, 
            eval_list = split(evalL)
            #print eval_list
            self.assertFloatEquals(float(strip(dsolvL)), round(float(strip(eval_list[2])), 3))


# AutogridTestCase


class BasicTestCase(Autogrid4TestCase):

    def test_testC_estat(self):
        """testC test estat case"""
        gpf_filename = 'testC.gpf'
        glg_filename = 'testC.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("C")
        self.compare_e_maps("C")


    def test_testC_dsolv(self):
        """testC test dsolv case"""
        gpf_filename = 'testC.gpf'
        glg_filename = 'testC.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("C")
        self.compare_d_maps("C")


    def test_testC_vdw(self):
        """testC test vdw case"""
        gpf_filename = 'testC.gpf'
        glg_filename = 'testC.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("C")
        self.compare_vdw_maps("C")


    def test_testA_estat(self):
        """testA test estat case"""
        gpf_filename = 'testA.gpf'
        glg_filename = 'testA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("A")
        self.compare_e_maps("A")


    def test_testA_dsolv(self):
        """testA test dsolv case"""
        gpf_filename = 'testA.gpf'
        glg_filename = 'testA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("A")
        self.compare_d_maps("A")


    def test_testA_vdw(self):
        """testA test vdw case"""
        gpf_filename = 'testA.gpf'
        glg_filename = 'testA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("A")
        self.compare_vdw_maps("A")


    def test_testN_estat(self):
        """testN test estat case"""
        gpf_filename = 'testN.gpf'
        glg_filename = 'testN.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("N")
        self.compare_e_maps("N")


    def test_testN_dsolv(self):
        """testN test dsolv case"""
        gpf_filename = 'testN.gpf'
        glg_filename = 'testN.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("N")
        self.compare_d_maps("N")


    def test_testN_vdw(self):
        """testN test vdw case"""
        gpf_filename = 'testN.gpf'
        glg_filename = 'testN.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("N")
        self.compare_vdw_maps("N")


    def test_testNA_estat(self):
        """testNA test estat case"""
        gpf_filename = 'testNA.gpf'
        glg_filename = 'testNA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("NA")
        self.compare_e_maps("NA")


    def test_testNA_dsolv(self):
        """testN test dsolv case"""
        gpf_filename = 'testNA.gpf'
        glg_filename = 'testNA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("NA")
        self.compare_d_maps("NA")


    def test_testNA_vdw(self):
        """testNA test vdw case"""
        gpf_filename = 'testNA.gpf'
        glg_filename = 'testNA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("NA")
        self.compare_vdw_maps("NA")


    def test_testOA_estat(self):
        """testOA test estat case"""
        gpf_filename = 'testOA.gpf'
        glg_filename = 'testOA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("OA")
        self.compare_e_maps("OA")


    def test_testOA_dsolv(self):
        """testO test dsolv case"""
        gpf_filename = 'testOA.gpf'
        glg_filename = 'testOA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("OA")
        self.compare_d_maps("OA")


    def test_testOA_vdw(self):
        """testOA test vdw case"""
        gpf_filename = 'testOA.gpf'
        glg_filename = 'testOA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("OA")
        self.compare_vdw_maps("OA")

    def test_testSA_estat(self):
        """testSA test estat case"""
        gpf_filename = 'testSA.gpf'
        glg_filename = 'testSA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("SA")
        self.compare_e_maps("SA")


    def test_testSA_dsolv(self):
        """testS test dsolv case"""
        gpf_filename = 'testSA.gpf'
        glg_filename = 'testSA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("SA")
        self.compare_d_maps("SA")


    def test_testSA_vdw(self):
        """testSA test vdw case"""
        gpf_filename = 'testSA.gpf'
        glg_filename = 'testSA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("SA")
        self.compare_vdw_maps("SA")


    def test_testHD_estat(self):
        """testHD test estat case"""
        gpf_filename = 'testHD.gpf'
        glg_filename = 'testHD.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("HD")
        self.compare_e_maps("HD")


    def test_testHD_dsolv(self):
        """testS test dsolv case"""
        gpf_filename = 'testHD.gpf'
        glg_filename = 'testHD.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("HD")
        self.compare_d_maps("HD")


    def test_testHD_vdw(self):
        """testHD test vdw case"""
        gpf_filename = 'testHD.gpf'
        glg_filename = 'testHD.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        #self.compare_maps("HD")
        self.compare_vdw_maps("HD")


    def test_testA(self):
        """testA test case"""
        gpf_filename = 'testA.gpf'
        glg_filename = 'testA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("A")


    def test_testN(self):
        """testN test case"""
        gpf_filename = 'testN.gpf'
        glg_filename = 'testN.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("N")


    def test_testNA(self):
        """testNA test case"""
        gpf_filename = 'testNA.gpf'
        glg_filename = 'testNA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("NA")


    def test_testS(self):
        """testS test case"""
        gpf_filename = 'testS.gpf'
        glg_filename = 'testS.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("S")


    def test_testSA(self):
        """testSA test case"""
        gpf_filename = 'testSA.gpf'
        glg_filename = 'testSA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("SA")


    def test_testOA(self):
        """testOA test case"""
        gpf_filename = 'testOA.gpf'
        glg_filename = 'testOA.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("OA")


    def test_testHD(self):
        """testHD test case"""
        gpf_filename = 'testHD.gpf'
        glg_filename = 'testHD.glg'
        # run autogrid4
        cmd_str = "%s -p %s -l %s" % \
                  (self.autogrid, gpf_filename, glg_filename)
        self.run_cmd(cmd_str)
        # compare resulting map with saved map
        self.compare_maps("HD")


# BasicTestCase

    def test_floatEqual(self):
        """ check precision: values within +/-.001 pass"""
        self.assertFloatEquals(0.001, 0.002)
        self.assertFloatEquals(0.011, 0.012)
        self.assertFloatEquals(0.021, 0.022)
        self.assertFloatEquals(0.031, 0.032)
        self.assertFloatEquals(0.041, 0.040)
        self.assertFloatEquals(0.051, 0.050)
        self.assertFloatEquals(0.061, 0.060)

 
if __name__ == '__main__':
    unittest.main()
