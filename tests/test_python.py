"""Pure-Python NETI against the certified grid: python3 tests/test_python.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
sys.path.insert(0, os.path.join(HERE, "..", "conformance"))
import neti_eckart
from conformance import run
raise SystemExit(0 if run(neti_eckart.kappa, tol=1e-12) else 1)
