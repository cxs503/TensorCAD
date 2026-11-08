# Triangle Intersection Verification Report

## Summary
Found a failing test in `tests/test_intersect.py::TestExactDifferential::test_thirds_lattice_matches_halfspace_clipping`.

## Issue
The float implementation `triangles_intersect()` disagrees with two exact routes (`triangles_intersect_exact` and `triangles_intersect_clip3d`) on 12 out of 600 random triangle pairs with third-based (non-dyadic) coordinates.

## Root Cause
Third-based coordinates (1/3, -1/3) are not exactly representable in binary floating point. The rounding errors accumulate through the intersection calculation and exceed the tolerance used in the point-in-triangle test.

## Failing Cases
First failing case:
- Triangle 1: ((-0.333..., 1.0, 1.0), (0.0, 1.0, 0.333...), (-0.333..., -0.333..., -1.0))
- Triangle 2: ((-1.0, 1.0, -1.0), (-1.0, 0.333..., 1.0), (0.0, 0.333..., 1.0))

The exact routes agree that these triangles intersect, but the float implementation returns False.

## Attempted Fixes
1. Increased tolerance from 1e-15 to 1e-12: no effect
2. Increased tolerance from 1e-12 to 1e-10: no effect

The issue is not just tolerance - it's a fundamental limitation of the floating-point algorithm when dealing with non-dyadic coordinates.

## Recommendation
This is a known limitation of floating-point arithmetic. The exact routes (rational arithmetic) should be used when precision is critical. The float implementation is suitable for mesh-level self-intersection detection where the geometry is well-conditioned.

## Test Status
- 38/39 tests pass in test_intersect.py
- 1 test fails: test_thirds_lattice_matches_halfspace_clipping
- All other test suites pass (1662 tests total)
