# Fillet kernel, independent verification

| check | measured | reference | rel. error | tolerance | passed |
|---|---|---|---|---|---|
| area r=0.05 | 7.997853982 | 7.997853982 | 0.00e+00 | 1e-14 | yes |
| perimeter r=0.05 | 11.91415927 | 11.91415927 | 1.49e-16 | 1e-14 | yes |
| I_uu r=0.05 | 10.65817801 | 10.65817801 | 1.67e-16 | 1e-13 | yes |
| I_vv r=0.05 | 2.664568127 | 2.664568127 | 0.00e+00 | 1e-13 | yes |
| area r=0.2 | 7.965663706 | 7.965663706 | 0.00e+00 | 1e-14 | yes |
| perimeter r=0.2 | 11.65663706 | 11.65663706 | 1.52e-16 | 1e-14 | yes |
| I_uu r=0.2 | 10.53534038 | 10.53534038 | 3.37e-16 | 1e-13 | yes |
| I_vv r=0.2 | 2.635281409 | 2.635281409 | 1.69e-16 | 1e-13 | yes |
| area r=0.4 | 7.862654825 | 7.862654825 | 1.13e-16 | 1e-14 | yes |
| perimeter r=0.4 | 11.31327412 | 11.31327412 | 0.00e+00 | 1e-14 | yes |
| I_uu r=0.4 | 10.16450254 | 10.16450254 | 5.24e-16 | 1e-13 | yes |
| I_vv r=0.4 | 2.551995263 | 2.551995263 | 1.74e-16 | 1e-13 | yes |
| area r=0.8 | 7.450619298 | 7.450619298 | 3.58e-16 | 1e-14 | yes |
| perimeter r=0.8 | 10.62654825 | 10.62654825 | 1.67e-16 | 1e-14 | yes |
| I_uu r=0.8 | 8.831924211 | 8.831924211 | 0.00e+00 | 1e-13 | yes |
| I_vv r=0.8 | 2.28372386 | 2.28372386 | 3.89e-16 | 1e-13 | yes |
| area r=0.999 | 7.14330861 | 7.14330861 | 2.49e-16 | 1e-14 | yes |
| perimeter r=0.999 | 10.28490212 | 10.28490212 | 3.45e-16 | 1e-14 | yes |
| I_uu r=0.999 | 7.931847653 | 7.931847653 | 1.01e-15 | 1e-13 | yes |
| I_vv r=0.999 | 2.119589758 | 2.119589758 | 4.19e-16 | 1e-13 | yes |
| brute-force area | 7.8627 | 7.862654825 | 5.75e-06 | 2e-03 | yes |
| brute-force centroid u | 2 | 2 | 2.96e-11 | 1e-04 | yes |
| brute-force I_uu | 10.16464747 | 10.16450254 | 1.43e-05 | 1e-04 | yes |
| brute-force I_vv | 2.552026512 | 2.551995263 | 1.22e-05 | 1e-04 | yes |
| L area r=0.1 | 2.991415927 | 2.991415927 | 0.00e+00 | 1e-14 | yes |
| L perimeter r=0.1 | 7.742477796 | 7.742477796 | 1.15e-16 | 1e-14 | yes |
| L area r=0.3 | 2.922743339 | 2.922743339 | 3.04e-16 | 1e-14 | yes |
| L perimeter r=0.3 | 7.227433388 | 7.227433388 | 2.46e-16 | 1e-14 | yes |
| L area r=0.5 | 2.785398163 | 2.785398163 | 0.00e+00 | 1e-14 | yes |
| L perimeter r=0.5 | 6.71238898 | 6.71238898 | 0.00e+00 | 1e-14 | yes |
| chordal error ratio at 8 chords (order 2) | 3.976929613 | 4 | 5.77e-03 | 1e-01 | yes |
| chordal error ratio at 16 chords (order 2) | 3.994220876 | 4 | 1.44e-03 | 1e-01 | yes |
| chordal error ratio at 32 chords (order 2) | 3.998554498 | 4 | 3.61e-04 | 1e-01 | yes |
| chordal error ratio at 64 chords (order 2) | 3.999638579 | 4 | 9.04e-05 | 1e-01 | yes |
| mesh volume error ratio per doubling (order 2) | 3.989901662 | 4 | 2.52e-03 | 5e-01 | yes |
| sweep volume r 0.2->0.6 | 23.55362818 | 23.55362818 | 0.00e+00 | 1e-14 | yes |
| sweep mesh volume r 0.2->0.6 | 23.55100522 | 23.55362818 | 1.11e-04 | 5e-04 | yes |
| sweep volume r 0.6->0.2 | 23.55362818 | 23.55362818 | 1.51e-16 | 1e-14 | yes |
| sweep mesh volume r 0.6->0.2 | 23.55100522 | 23.55362818 | 1.11e-04 | 5e-04 | yes |
| sweep volume r 0.1->0.999 | 23.04896964 | 23.04896964 | 0.00e+00 | 1e-14 | yes |
| sweep mesh volume r 0.1->0.999 | 23.04338121 | 23.04896964 | 2.42e-04 | 5e-04 | yes |
| sweep volume r 0.4->0.4 | 23.58796447 | 23.58796447 | 0.00e+00 | 1e-14 | yes |
| sweep mesh volume r 0.4->0.4 | 23.58554328 | 23.58796447 | 1.03e-04 | 5e-04 | yes |
| max radius of a 2x2 square | 1 | 1 | 2.22e-16 | 1e-12 | yes |
| circle-limit area | 3.141592654 | 3.141592654 | 2.83e-16 | 1e-14 | yes |
| circle-limit perimeter | 6.283185307 | 6.283185307 | 2.83e-16 | 1e-14 | yes |
| over-limit radius refused with the limit named | 1 | 1 | 0.00e+00 | 0e+00 | yes |
| sliver_4x0.2_r0.0999 sampled area below exact (chords cut the arc) | 0.7912320138 | 0.7914330861 | 2.54e-04 | 5e-03 | yes |
| rectangle_r0.999 sampled area below exact (chords cut the arc) | 7.123201383 | 7.14330861 | 2.81e-03 | 5e-03 | yes |
| l_section_r0.5 sampled area below exact (chords cut the arc) | 2.780361288 | 2.785398163 | 1.81e-03 | 5e-03 | yes |
| notched_section_r0.29 sampled area below exact (chords cut the arc) | 2.326113537 | 2.327807942 | 7.28e-04 | 5e-03 | yes |
| clockwise input area equals counter-clockwise | 7.862654825 | 7.862654825 | 0.00e+00 | 1e-14 | yes |
| area error at a 1e6 translation (1e6/4 = 2.5e5 x the body) | 2.842781451e-11 | 0 | 2.84e-11 | 1e-09 | yes |
