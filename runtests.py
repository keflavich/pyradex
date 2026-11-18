#!/usr/bin/env python
"""
Simple test runner using pytest directly.
This replaces the old deprecated astropy test runner.
"""
import sys
import pytest

if __name__ == '__main__':
    # Run pytest on the tests directory
    sys.exit(pytest.main(['pyradex/tests/'] + sys.argv[1:]))
