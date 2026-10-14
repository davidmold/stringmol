# Running Stringmol in Docker

Branch `perf`: upstream commit 8e568df plus performance fixes that keep output
byte-identical (see git log). The image pins the toolchain (CUDA 13.2 devel,
g++ 13.3) by digest so reference and future GPU builds share one compiler.

```
./sm build
./sm start ref_v2.0.3 config/spatial_stringmol_v2.0.3.conf
./sm status
./sm logs ref_v2.0.3
./sm stop ref_v2.0.3
```

Each run lives in `runs/NAME/` (run.conf, splist*.dat, out1_*.conf, popdy001.dat).
Progress lines go to the container log (`./sm logs`).

Notes
- CPU-only, single-threaded per run; no GPU reservation. aida has 6 physical
  cores, so keep concurrent runs to about 6.
- A full 1M-step run writes ~10 GB (a full state dump every 100 steps).
- The step loop is hard-coded to 1,000,000 steps in `src/smspatial.cpp`; NSTEPS
  in the config is ignored. Step 90,001 prints a block of RNG debug lines.
- Probes: `./sm bind probe/pairs.txt` and `./sm react probe/pairs3.txt` (binding table; full reactions with mutation off).
- RANDSEED is honoured from commit "honour RANDSEED" onward. Upstream never
  seeded the RNG in smspatial, so every earlier run (including ref_v2.0.3)
  used the Mersenne Twister default seed; reproduce those with RANDSEED 4357.
