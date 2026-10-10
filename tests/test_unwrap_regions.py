# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Water-mask regions (PRD R-U1; plan T36.1)."""

from __future__ import annotations

import numpy as np
import pytest

from venti.spatial.resample import downsample_array
from venti.unwrap.regions import downsample_labels, largest_region, segment_regions


def archipelago():
    """Mainland, an island, a peninsula on a 4-px neck, a sliver, an islet, a hole."""
    land = np.zeros((200, 300), bool)
    land[:, :120] = True  # mainland
    land[40:90, 160:230] = True  # island A, 40-px channel
    land[120:124, 120:150] = True  # neck, 4 px wide
    land[100:160, 150:210] = True  # peninsula B on the neck
    land[10:13, 250:253] = True  # 9-px islet
    land[180:184, 200:290] = True  # a 4-px-wide sliver (360 px), gone after erosion
    valid = np.ones(land.shape, bool)
    valid[60:70, 20:40] = False  # a no-data hole inside the mainland
    return land, valid


def test_regions_follow_the_water_mask():
    land, valid = archipelago()
    labels = segment_regions(land, valid)
    assert labels.dtype == np.int32
    # mainland + peninsula (joined by coherent land) and the island; the islet
    # (< 20 px) and the sliver (no seed survives three erosions) are dropped
    assert labels.max() == 2
    assert set(np.unique(labels)) == {0, 1, 2}
    main = largest_region(labels)
    assert np.all(labels[:, :100][valid[:, :100]] == main)
    assert labels[11, 251] == 0  # islet
    assert np.all(labels[180:184, 200:290] == 0)  # sliver
    assert len(np.unique(labels[50:80, 170:220])) == 1  # one label on the island
    island = int(labels[60, 190])
    assert island != main
    assert labels[140, 180] == main
    assert labels[122, 135] == main
    # water and no-data are 0
    assert labels[100, 140] == 0
    assert np.all(labels[60:70, 20:40] == 0)


def test_erosion_decides_what_is_too_thin_and_parameters():
    land = np.zeros((100, 200), bool)
    land[:, :80] = True
    land[40:50, 120:180] = True  # a 10-px-wide strip 40 px off the coast
    assert segment_regions(land).max() == 2  # survives 3 erosions
    assert segment_regions(land, erosion_iterations=6).max() == 1  # not 6
    assert segment_regions(land, erosion_iterations=0).max() == 2
    land[40:50, 80:120] = True  # join it to the mainland by a 10-px neck
    assert segment_regions(land).max() == 1
    with pytest.raises(ValueError, match="2-D"):
        segment_regions(land[0])
    with pytest.raises(ValueError, match="valid_mask"):
        segment_regions(land, np.ones((3, 3), bool))
    with pytest.raises(ValueError, match="min_region_area"):
        segment_regions(land, min_region_area=0)
    assert largest_region(np.zeros((4, 4), np.int32)) == 0


def test_downsample_labels_matches_the_fit_grid():
    land, valid = archipelago()
    labels = segment_regions(land, valid)
    for factor in (1, 3, 6, 7):
        small = downsample_labels(labels, factor)
        ref = downsample_array(labels.astype(np.float32), factor)
        assert small.shape == ref.shape
        if factor == 1:
            assert small is labels
    small = downsample_labels(labels, 6)
    # the island keeps its label and the water stays 0
    assert small[60 // 6, 190 // 6] == labels[60, 190]
    assert small[100 // 6, 140 // 6] == 0
    with pytest.raises(ValueError, match="factor"):
        downsample_labels(labels, 0)
