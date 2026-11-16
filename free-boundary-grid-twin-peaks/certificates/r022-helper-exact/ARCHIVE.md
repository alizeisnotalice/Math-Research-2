# r022 corrected helper exact-DP certificate

Source snapshot: `a022-helper-exact-dp-r022`.  This package retains only the
corrected n=3/n=4 complete results (`n3.json`, `n4.json`) and the final C++
source plus independent n=3 direct audit.  The half-cover outputs are
intentionally absent: they are known buggy artifacts and are not evidence.

The recorded exact totals are n=3: 12 cover records, weight 129024, sum
g/h=51072; n=4: 72 cover records, weight 1194718147840, sum
g/h=975017025920/3.  Build outside this archive if desired:
`c++ -O3 -std=c++20 helper_dp.cpp -o helper_dp`, then run `./helper_dp 3`
or `./helper_dp 4`; use `python3 -B audit_n3_direct.py` for the independent
n=3 audit.  The executable is deliberately not archived.
Final corrected C++ SHA-256: `be0ae1934efa8622d1281fbb27a833c94bc99b00dfc3f82df78024efe2d1375e`.
