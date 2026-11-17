# R19 independent native-Metal orbit verification receipt

Input: `orbits-gpu.json`, SHA-256
`feba42b40fa4185cba8062ddbf662221ab889418028fc8775bc1681fb0a3eab7`.
Source: `R19.swift`, SHA-256
`4dc2cf8030d399eea7b10e5d365452fb164e24095e5365a17275589d0fbdd787`.
Result: `result.json`, whose embedded source/input hashes bind those exact
files and whose weighted count is `232081500`.

Remote directory was
`/Users/zhengzhihao/CodexRemoteMath/a28-selector/r019/metal-independent` on
the MacBook.  The actual foreground command was equivalent to:

```sh
ih=$(shasum -a 256 orbits-gpu.json | awk '{print $1}')
sh=$(shasum -a 256 R19.swift | awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift-module-cache" \
  -framework Foundation -framework Metal R19.swift -o r19
./r19 orbits-gpu.json result.json "$ih" "$sh"
```

The same MacBook environment used for the surrounding R16--R18 Metal runs was
Apple Swift 6.4.2 on macOS 27.3 arm64; this R19 task did not separately capture
an environment file, so no stronger per-run version timestamp is claimed.
Start/end timestamps and elapsed host time were not instrumented for R19 and
are intentionally not reconstructed.  The SSH foreground command returned
successfully, and a subsequent process query showed no owned `r19` or `swiftc`
process, so its Metal context and buffers had exited.

Every one of 6 orbit representatives and 32,768 quadratic phases was sent to
the native 64-point int32 Metal FWHT.  CPU int64 FWHT compared all signed
coefficients without retaining large output arrays:
`64 * 6 * 32768 = 12,582,912` coefficient comparisons.  The Metal output
buffer was at most `2048 * 64 * 4 = 524,288` bytes (under the 2 MiB limit).
The per-orbit fiber counts were `2604, 784, 384, 448, 0, 0`, agreeing with the
separate CPU result; weighting by the supplied orbit sizes gives `232081500`.

`R19.swift` prints a report whose `status` is authoritative.  Its final Swift
line currently does not call `exit(1)` if that status were `failed`; therefore
the successful shell exit is supplemental evidence, while reproducible
acceptance must check `result.json.status == "passed"`.  No source change or
rerun was made merely to alter that behavior.
