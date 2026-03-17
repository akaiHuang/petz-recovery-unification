"""Centralized type aliases for tau-chrono."""

from __future__ import annotations

from typing import List

import numpy as np
from numpy.typing import NDArray

DensityMatrix = NDArray[np.complexfloating]
KrausList = List[NDArray[np.complexfloating]]
SuperOp = NDArray[np.complexfloating]
