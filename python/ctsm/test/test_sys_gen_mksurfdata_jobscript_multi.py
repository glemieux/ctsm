#!/usr/bin/env python3

"""System tests for gen_mksurfdata_jobscript_multi

"""

import os

import unittest
import tempfile
import shutil
import sys

from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_jobscript_multi import main
from ctsm import unit_testing

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysGenMkSurfJSMulti(unittest.TestCase):
    """System tests for gen_mksurfdata_jobscript_multi"""

    def setUp(self):
        """Setp temporary directory to make the files in"""
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)
        self.outfile = "jobscript.sh"
        bld_path = os.path.join(
            path_to_ctsm_root(), "python", "ctsm", "test", "testinputs", "mksurfdata_esmf_bld"
        )
        sys.argv = [
            "gen_mksurfdata_jobscript_multi",
            "--bld-path",
            bld_path,
            "--jobscript-file",
            self.outfile,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def createJS(self, nodes, tasks_per_node, scenario, option_list=[]):
        """
        Create a JobScript by sending a list of options in
        """
        if len(option_list) > 1:
            sys.argv.extend(option_list)
        sys.argv.extend(
            [
                "--number-of-nodes",
                nodes,
                "--tasks-per-node",
                tasks_per_node,
                "--scenario",
                scenario,
            ]
        )
        print(sys.argv)
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output jobscript file should exist")

    def test_simple_jobscript_multi(self):
        """
        Test that a standard simple namelist works
        """
        # pylint: disable=no-self-use
        self.createJS(nodes="4", tasks_per_node="12", scenario="standard")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
