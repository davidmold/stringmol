# Stringmol, fast

A maintained fork of [Stringmol](https://github.com/franticspider/stringmol), the automata
chemistry by Simon Hickinbotham, Susan Stepney and colleagues at the University of York.
In Stringmol, molecules are strings of opcodes that bind by sequence alignment and then
run as programs that copy, cut and edit each other. The papers listed below describe the
chemistry; `README-upstream.md` is the original README.

This fork targets the **spatial** simulator (`smspatial`). The chemistry is unchanged:
every change below was checked to give **byte-identical results** to the original for the
same seed. What changes is speed, output size, a few real bugs, and some tools for looking
at runs.

## Why use this fork

| | upstream (8e568df) | this fork |
|---|---|---|
| speed, early run (~5,000 molecules) | ~4 steps/s at step 20,000 | ~240 steps/s |
| speed, full 125x100 grid (~11,000 molecules) | under 1 step/s | ~120 steps/s |
| disk per 1M-step run | ~100 GB | a few hundred MB |
| `RANDSEED` in the config | ignored (always seed 4357) | honoured |
| restart from a saved state | silently reseeds randomly | reproduces the run |
| extinct runs | continue to step 1M on an empty grid | stop |

Timings are single-threaded on a desktop CPU and depend on what the population is doing.

The fork starts from upstream commit `8e568df`. The three upstream commits after it (a
memory-leak cleanup) make both simulators segfault, so they are not included.

## Build and run

Linux, g++ and make:

```
cd src && mkdir -p ../release && make all
mkdir -p ../runs/demo && cd ../runs/demo
cp ../../config/spatial_stringmol_v2.0.3.conf run.conf
../../src/smspatial x run.conf
```

Or with Docker (see `DOCKER.md`): `./sm build`, then
`./sm start demo config/spatial_stringmol_v2.0.3.conf`, `./sm status`, `./sm viz demo`.

Set `RANDSEED` in the config to get a different run. Runs made with upstream builds are
reproduced with `RANDSEED 4357`.

## Output

Written into the directory the simulator runs in:

| file | contents |
|---|---|
| `popdy001.dat` | `step,species,count` for every living species, every 100 steps (same format as upstream) |
| `species.tsv` | each species once, when it first appears: id, birth step, parents, sequence |
| `grid.snap` | the grid every `SNAPEVERY` steps (default 500): per frame a 16-byte header (step, width, height, byte count) then zlib-compressed `uint32` species ids, row-major, 0 for empty, bit 31 set for bound molecules |
| `out1_N.conf`, `RNGstate_N.txt` | full checkpoints every `CKPTEVERY` steps (default 50,000); the newest `KEEPCKPT` (default 3) are kept |
| `splist_final.dat` | the full upstream-format species list, once, at the end |

`LEGACYOUT 1` in the config brings back the upstream output (species list and full state
every 100 steps) if you need it. Expect around 100 GB per million steps.

To resume a run, copy a checkpoint's `out1_N.conf` and `RNGstate_N.txt` into a new
directory and run `smspatial x out1_N.conf` there. The population dynamics continue exactly;
species that went extinct and later re-appear get new id numbers.

## Tools

- **Viewer** (`viz/build.py RUN_DIR`): writes a self-contained `viewer.html` with the grid
  animated over time, the most abundant lineages coloured, and a population timeline you
  can scrub. Works on both slim and legacy output.
- **Binding probe** (`probe/probe.cpp`): binding probability and copier/copied roles for
  pairs of molecules, using the simulator's own binding code.
- **Reaction probe** (`probe/react.cpp`): runs a reaction between two molecules to
  completion, mutation off, and reports what it produced.

Both probes are built into the Docker image (`./sm bind FILE`, `./sm react FILE`); the
pairs-file format is in `sm`.

## What changed, and how it was checked

Performance, all output-preserving:
- cached opcode lookup in the alignment scoring function;
- the agent lists keep an ordered mirror (a Fenwick tree over append-only slots) so
  append, removal and picking the k-th molecule are O(log n) while list order, and so the
  random-number stream, is unchanged;
- species lookup through a hash index instead of a scan of every species ever seen;
- the population writer makes one pass instead of one pass per species;
- a Smith-Waterman alignment whose result was discarded is no longer computed.

Bugs fixed:
- `smspatial` never seeded the random number generator, so `RANDSEED` had no effect;
- the path to `RNGFILE` was built wrongly for configs in the current directory, so restarts
  could not load the saved generator state and reseeded from `/dev/urandom`;
- runs carried on to step 1M after the population died out.

Each change was checked against the previous build by comparing `popdy001.dat` and the
species list byte for byte over the first 20,000+ steps, and the list-mirror changes also
over a resumed full-grid run (steps 400,000-410,000). Setting `STRINGMOL_VERIFY=1` makes the
simulator check the list mirrors against the linked lists on every access and report any
drift (slow; for testing only).

Known quirks left as they are, because they are part of the chemistry: the binding
probability normalises alignment score by length, so very short fragments bind easily
(the source calls this a "BRUTAL HACK"); the step loop is fixed at 1,000,000 steps and
`NSTEPS` is ignored; step 90,001 prints a block of RNG debug lines.

## Licence

GPL, as upstream. Note that upstream's `LICENSE` file is GPL version 2 while the source
file headers say "GPL version 3 or (at your option) any later version"; both are kept as
they are. Copyright in the original code remains with its authors.

## Papers

- Hickinbotham et al., *Specification of the Stringmol chemical programming language*,
  University of York tech report YCS-2010-458 (`smtr0.2.pdf` in this repo).
- Hickinbotham et al., *Diversity from a monoculture: effects of mutation-on-copy in a
  string-based artificial chemistry*, ALife XII, 2010.
- Clark et al., *Degeneracy enriches artificial chemistry binding systems*, ECAL 2011.
- Hickinbotham et al., *Maximizing the adjacent possible in automata chemistries*,
  Artificial Life 22(1), 2016.
- Hickinbotham & Stepney, *Conservation of matter increases evolutionary activity*,
  ECAL 2015.
- Hickinbotham, Stepney & Hogeweg, *Nothing in evolution makes sense except in the light of
  parasites*, 2021.
