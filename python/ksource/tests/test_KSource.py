#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

import numpy as np
from ksource import KSource

@pytest.mark.parametrize("number", [1])
def test_test(number):
	print("Python tests still not implemented.")

if __name__ == "__main__":
    # --durations=10  <- May be used to show potentially slow tests
    pytest.main(args=[__file__, "--doctest-modules", "-v", "--durations=15"])

