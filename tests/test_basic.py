# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Basic smoke tests that don't depend on other venti modules.

Behaviour of the unwrap corrector is tested in test_unwrap_corrections.py.
"""

import numpy as np
import pytest


def test_numpy_available():
    """Test that numpy is available."""
    arr = np.array([1, 2, 3])
    assert arr.sum() == 6


def test_can_import_unwrap():
    """Test that we can import the unwrap module without triggering pyproj."""
    # This should not trigger models import
    from venti.unwrap import UnwrapCorrector

    assert UnwrapCorrector is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
