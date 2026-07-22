# Environment used to build and debug MPI + adjoint (2026-07-20/21)

This documents the exact toolchain that was used to find and fix the MPI/adjoint
bugs in `debug_mpi/` (heap corruption in `com_dof`, silently-dropped adjoint of
`mpi_sum_r`/`mpi_sum_i`, cost double-counted under MPI). Several of the bugs
were tooling-version-sensitive (see "Known fragile points" below) rather than
pure DassFlow logic bugs, so pinning/recording these versions matters for
reproducing a working build later, on another machine, or after a system
upgrade.

## OS / system toolchain

| Component | Version |
|---|---|
| OS | Ubuntu 22.04.5 LTS (jammy) |
| Kernel | 6.8.0-134-generic |
| gcc | 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3) |
| gfortran | 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3) |
| OpenMPI | 4.1.2-2ubuntu1 (apt: `libopenmpi-dev`, `libopenmpi3`, `openmpi-bin`, `openmpi-common`) |
| Java (for Tapenade) | OpenJDK 11.0.31 |
| Perl (for `finish_to_gen_adjoint.pl` / `extract_df_sum_adjoint.pl`) | 5.34.0 |

Install the apt packages with:
```
sudo apt install build-essential gfortran libopenmpi-dev openmpi-bin openjdk-11-jdk perl
```

## Vendored / bundled libraries (in `code/libs/`, built from source by `make lib`)

| Library | Version | Notes |
|---|---|---|
| SCOTCH | 5.1.12 (esmumps variant) | `code/libs/scotch_5.1.12_esmumps/` |
| MUMPS | 4.10.0 | `code/libs/MUMPS_4.10.0/` |
| m1qn3 | 3.3 | `code/libs/m1qn3-3.3-distrib/` |
| AGMG | 3.2.0-aca | `code/libs/AGMG_3.2.0-aca/` |
| Tapenade ADFirstAidKit | vendored 2026-07-20 from Tapenade 3.16 | `code/libs/adFirstAidKit/` -- **new**, added for the MPI-adjoint fix, see below |

These build against the system gcc/gfortran/OpenMPI above -- if you change any
of those, rebuild them too (`make cleanlib; make lib`).

## Tapenade (source-to-source AD)

```
Tapenade 3.16 (develop) - 16 Jul 2026 13:00 - Java 11.0.31 Linux
Revision: e59864cab441d4175df75383b3ff58c3dcd26df9
Tag:      3.16-v2-1404-ge59864cab
```

This is a `develop` snapshot, not a tagged release -- there is no guarantee an
older or newer Tapenade build differentiates `df_sum_r`/`df_sum_i` the same
way (see fragile point below). If you need to re-fetch Tapenade, prefer this
exact revision, or re-verify the MPI-adjoint fix (`code/src/adjoint/extract_df_sum_adjoint.pl`,
`code/src/common/m_mpi.f90`) against whatever version you land on.

## Python (conda environment `dassflow2`)

```
$ LOAD_LIB
$ conda activate dassflow2
$ python3 --version
Python 3.9.12
```

`requirements.txt` (next to this file) is the exact `pip freeze` of this
environment. The ones that actually matter for building/running DassFlow2D:

| Package | Version | Why it matters |
|---|---|---|
| f90wrap | **0.2.16** | See fragile point below -- must match what the `.so` was built against. |
| numpy | 2.0.2 | |
| scipy | 1.13.1 | |
| mpi4py | 4.1.2 | Only used to trigger `MPI_Init`/`MPI_Finalize` from Python; must be built against the same OpenMPI 4.1.2 as everything else. |
| pandas, matplotlib, h5py, pyvista, vtk | see requirements.txt | Post-processing / plotting, not build-critical. |
| setuptools | whatever pip resolves | See fragile point below (distutils shadowing). |

Note: `dassflow2d/environment.yml` (repo root, pre-existing) pins
`python=3.8` and `f90wrap==0.2.7` -- that is **not** what this working build
actually uses (Python 3.9.12 / f90wrap 0.2.16). Treat `environment.yml` as
stale/aspirational until someone reconciles it with `requirements.txt`.

## Known fragile points (found the hard way during this debugging session)

1. **f90wrap version is an ABI contract with the compiled `.so`, not just a
   Python API.** Installing a different f90wrap (e.g. a bare `pip install
   f90wrap` pulling latest, 0.3.0) into the *same* interpreter search path
   (even via `~/.local`, which Python prepends ahead of the conda env's own
   site-packages) silently shadows the working 0.2.16 and produces bare,
   message-less `ValueError`s deep in `arraydata`/derived-type array access,
   with no indication it's a version mismatch. If you ever see mysterious
   `ValueError()` with an empty message from a `wrapping/m_*.py` getter,
   check `python -c "import f90wrap; print(f90wrap.__file__, f90wrap.__version__)"`
   resolves to the conda env, not `~/.local`.

2. **Recent `setuptools` (>=~60) hijacks the stdlib `distutils`**, which
   breaks `numpy.f2py`'s build backend when `f90wrap`'s `f2py-f90wrap` CLI
   invokes it (`ModuleNotFoundError: No module named
   'distutils.msvccompiler'`). Workaround: `export
   SETUPTOOLS_USE_DISTUTILS=stdlib` before `make install`/`make wrap`.

3. **Tapenade's own Fortran parser needs its own `-I` for `mpif.h`,
   separately from the compiler's.** `df_sum_r`/`df_sum_i` (`m_mpi.f90`)
   deliberately show Tapenade the real `include 'mpif.h'` +
   `MPI_ALLREDUCE` call (no differentiable stub) so it auto-differentiates
   correctly via its built-in `MPI_Allreduce_fwd/bwd/d` support. If Tapenade
   can't resolve `mpif.h` itself (it looks relative to `tap/`, not via the
   compiler's include path), it corrupts its internal AST and dies with an
   unrelated-looking `java.lang.StackOverflowError` much later, not a
   parse error. The Makefile computes this via `mpicc -show` and passes it
   as `-I` to `tapenade`, plus `-java "-Xss64m"` for extra JVM stack
   headroom -- if `mpicc`/OpenMPI's include layout ever changes, re-check
   `TAPENADE_MPI_I` in the Makefile still resolves to a real path.

4. **`finish_to_gen_adjoint.pl` deletes every generated `*mpi_back*`/
   `*mpi_diff*` file wholesale** (to discard Tapenade's own wrong
   differentiation of `com_dof`/`com_var_r`, replaced by the hand-written
   versions in `src/adjoint/m_mpi_back.f90`). Since `df_sum_r`/`df_sum_i`'s
   *correct* auto-generated adjoint also lives in that doomed file,
   `code/src/adjoint/extract_df_sum_adjoint.pl` runs first and rescues just
   those two subroutines into `df_sum_reduce_back.f90`/`_diff.f90`. If
   Tapenade's output shape for `DF_SUM_R_BACK`/`DF_SUM_R_DIFF` ever changes
   (different Tapenade version), this extraction script's pattern match may
   need adjusting -- it prints a `WARNING ... adjust this script` if it
   can't find what it's looking for, don't ignore that warning.

5. **`adMPI.c` (vendored in `code/libs/adFirstAidKit/`) calls a handful of
   push/pop primitives as plain C functions** (`pushCharacterArray`,
   `popCharacterArray`, `pushInteger4`, `popInteger4`), while DassFlow's
   existing stack (`src/adjoint/adBuffer.f`) only exposes the Fortran-mangled
   equivalents. `code/libs/adFirstAidKit/adCharacterArrayShim.c` bridges the
   four that are actually used -- re-derive this list with
   `grep -oE '\b(push|pop|look)[A-Za-z0-9]*\(' code/libs/adFirstAidKit/adMPI.c | sort -u`
   if `adMPI.c` is ever updated to a newer Tapenade version, since it may
   call additional ones.

6. The shell profile sets `export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6`
   globally (in `start_system.sh`, loaded by the `LOAD_LIB` alias). Not
   something this fix depends on, but worth knowing it's there if you ever
   see unexplained C++ ABI behavior.
