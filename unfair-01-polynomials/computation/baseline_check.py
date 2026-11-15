#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import sys
import unittest
import numpy
import scipy
import sympy
from test_exact import ExactChecks
from algebraic_verify import verify

suite=unittest.defaultTestLoader.loadTestsFromTestCase(ExactChecks)
result=unittest.TextTestRunner(verbosity=2).run(suite)
data={'python':sys.version,'executable':sys.executable,'sympy':sympy.__version__,
      'numpy':numpy.__version__,'scipy':scipy.__version__,
      'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
      'exact_example':verify({'minpoly':[0,1],'interval':[-1,1],
                              'A':[[1],[1]],'B':[[1],[0],[1]]}),
      'script_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.glob('*.py'))}}
print(json.dumps(data,indent=2))
sys.exit(0 if result.wasSuccessful() else 1)
