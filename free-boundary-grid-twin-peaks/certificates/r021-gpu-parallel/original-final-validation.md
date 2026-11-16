# Final validation receipt

The final source SHA-256 is
`6aa1f829dbe97c141781c7440c0f99ef62fa966cbe2c665fab7432a48b2ec9f8`;
the final Swift binary SHA-256 is
`376b39276c96ed9e783a7edfc780fd40026cc94de42180d92b13f01b0e86fb6f`.

After the 24-bit entry guard was added, a real n=3, P=17, tau=10, L=16
trace again agreed exactly with r020 in every mathematical output field and
in each GPU m2/m3 group.  The independent terminal Python DP agreed too.

An intentional illegal call with n=2 and P=16,777,217 exited with status 133
at the Swift `precondition`, before output-file creation and before `run()`.
This validates the requested rejection ordering; it is not a GPU failure.
No `inverse_parallel` process remained after these checks.
