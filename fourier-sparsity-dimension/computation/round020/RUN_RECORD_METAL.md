# R20 two-host native Metal verification

The root actually compiled and ran final R20Metal.swift on both machines.
Final source SHA256:40d770e83346d922bb0c592eb799d0cc2f8cf5c5248e1e8a50252d7cd86a5abf.
Final input SHA256:5e7303916cc3029c7088955e5bed677071878cf717ed4578f078c46d02641789.
Swift6.4.2 (swiftlang-6.4.2.1.9), local macOS27.2 and remote macOS27.3.

Both runs first passed smoke at N16, then all120 pairs/240 phases, N16..4096.
Each run compared261888 signed int32 GPU coefficients to independent Swift
int64 CPU FWHT and checked magnitude q. Both report zero mismatches and
zero wrong magnitudes. Local device Apple M6; MacBook device Apple M5 Pro.
Reported elapsed times are0.022400021553039551 s and0.10607504844665527 s.
These timers include CPU transforms, GPU dispatch/wait/readback and comparison;
they start after Metal library compilation and input parsing, not end-to-end
startup. No sustained utilization or performance superiority is claimed.

Implementation: separate compute encoder for each butterfly stage, ping-pong
shared buffers with default tracked hazards. No oversized threadgroup is
required. Each f/g batch contains24 phases; the largest single buffer is
24*4096*4=393216 bytes. All operations are exact integers; zero tolerance.

Actual compilation: swiftc -O -module-cache-path CACHE -framework Foundation
-framework Metal R20Metal.swift -o R20Metal. Local CACHE was
/Volumes/data/CodexTemp/cache/swift-module-cache; remote CACHE was
/Users/zhengzhihao/CodexRemoteMath/cache/swift-module-cache.
Actual invocation in each dedicated directory:

```sh
./R20Metal smoke cases-gpu-final.json smoke.json
./R20Metal full cases-gpu-final.json result.json
```

The root's final local source fixed a Swift raw-string closing delimiter
that caused the first compile to fail. Before execution it also required
ceil(N/64) input word length and a real CryptoKit SHA256 input binding.
Only corrected, actually rerun results are archived. Intermittent SSH banner
timeouts occurred before final dispatch; a read-only remote check confirmed
no existing final GPU job, then the root uploaded the verified source and
ran smoke/full through one foreground SSH command. Compilation and both
programs exited0 on each host. No background GPU loop or MATLAB session was
started. The remote result was retrieved and its input/source bindings checked.
