# r023 fixed-state root completion

Retains final D4 source, subset-reference source, n=3..7 exact outputs,
small reference outputs, provenance and resources. `prepare.py` is preserved
as historical provenance even though it has original paths; compile and run
the relative sources directly, e.g. `c++ -O3 -std=c++20 fixed_state_d4.cpp
-o fixed_state_d4 && ./fixed_state_d4 5`. No binary is archived.
