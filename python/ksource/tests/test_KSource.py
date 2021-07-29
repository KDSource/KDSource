#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

import numpy as np
from ksource import KSource

@pytest.mark.parametrize("number", [1])
def test_test(number):
	assert number == 1

if __name__ == "__main__":
    # --durations=10  <- May be used to show potentially slow tests
    pytest.main(args=[__file__, "--doctest-modules", "-v", "--durations=15"])

