"""Pure-Python ELIAH against the certified grid: python3 tests/test_python.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
sys.path.insert(0, os.path.join(HERE, "..", "conformance"))
import eliah
from conformance import run
raise SystemExit(0 if run(eliah.kappa, tol=1e-12) else 1)
