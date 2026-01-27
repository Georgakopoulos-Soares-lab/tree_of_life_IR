import os
import sys

# Ensure repository root is on sys.path for package-less imports
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Also add the nonbdna_pipeline folder so absolute sibling imports inside that
# package (e.g., `from stream_and_merge_bucket import StreamAndMerge`) resolve.
NONBDNA = os.path.join(ROOT, "nonbdna_pipeline")
if NONBDNA not in sys.path:
    sys.path.insert(0, NONBDNA)
