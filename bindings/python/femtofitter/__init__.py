#
# femtofitter/__init__.py
#

import os
from pathlib import Path

LIBFEMTOFITTER_PATH = os.environ.get("LIBFEMTOFITTER")

if not LIBFEMTOFITTER_PATH:
    libfem = Path(__file__).parent
else:
    libfem = Path(LIBFEMTOFITTER_PATH)




from ROOT import gSystem
