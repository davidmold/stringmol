# Running Stringmol in Docker

The image pins the toolchain (CUDA 13.2 devel,
g++ 13.3) by digest so reference and future GPU builds share one compiler.

```
./sm build
./sm start ref_v2.0.3 config/spatial_stringmol_v2.0.3.conf
./sm status
./sm logs ref_v2.0.3
./sm stop ref_v2.0.3
```

Each run lives in `runs/NAME/` (see README.md for the output files). Progress lines go
to the container log (`./sm logs`). `./sm resume NEW OLD STEP` continues run OLD from
its checkpoint at STEP.

Notes
- CPU-only, single-threaded per run; no GPU reservation. Run about one
  simulation per physical core.
- Output (slim, default): popdy001.dat every 100 steps; species.tsv (each species once:
  id, birth step, parents, sequence); grid.snap (compressed species-id grid every
  SNAPEVERY=500 steps, ~2 KB/frame); out1_N.conf + RNGstate_N.txt checkpoints every
  CKPTEVERY=50000 steps, newest KEEPCKPT=3 kept; splist_final.dat at the end.
  A 1M-step run is a few hundred MB. LEGACYOUT 1 in the config restores the upstream
  every-100-steps dumps (~100 GB per run).
- The step loop is hard-coded to 1,000,000 steps in `src/smspatial.cpp`; NSTEPS
  in the config is ignored. Step 90,001 prints a block of RNG debug lines.
- Probes: `./sm bind probe/pairs.txt` and `./sm react probe/pairs3.txt` (binding table; full reactions with mutation off).
- RANDSEED is honoured from commit "honour RANDSEED" onward. Upstream never
  seeded the RNG in smspatial, so every earlier run (including ref_v2.0.3)
  used the Mersenne Twister default seed; reproduce those with RANDSEED 4357.
- Viewer: `./sm viz NAME` (or `python3 viz/build.py RUN_DIR`) writes a self-contained
  `viewer.html` into the run folder: the grid animated through time, top lineages
  coloured, and a scrubbable population timeline. Rebuild it to pick up new steps.
