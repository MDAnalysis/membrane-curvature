# MembraneCurvature benchmarks

Performance benchmarks for MembraneCurvature using [Airspeed Velocity] (ASV).

ASV creates isolated environments with the [uv plugin] (See [environment_type: uv] in [asv.conf.json]).

Each benchmark run uses Python **3.11** and **3.14**, the minimum and latest versions declared in [pyproject.toml].

## Prerequisites

- [uv] on your `PATH` (`uv --version`)

Install the required Python versions in your `uv` environment:

```bash
uv python install 3.11 3.14
```

From the repository root, sync the benchmark dependencies:

```bash
uv sync --extra benchmarks
```

If this is the first time you are running the benchmarks, register the machine with:

```bash
cd benchmarks
uv run asv machine --yes
```

## Usage

> [!NOTE]
> Run all `asv` commands from the `benchmarks/` directory.

### Benchmark the current commit

Benchmarks `HEAD` on the current checkout against the `main` branch. Results are saved, but the HTML only generates data if that commit is on a branch listed in `asv.conf.json` in the `main` branch.

```bash
cd benchmarks
uv run asv run HEAD --steps 1
uv run asv publish
uv run asv preview --browser
```

Timings are saved under `results/`, but summary graphs in the published site only include branches listed in [asv.conf.json] (currently `main`).

> [!WARNING]
> Do not open `html/index.html` via `file://`. Use `asv preview` to view the results in a local HTTP server.

### History on `main` (recommended for the benchmark grid)

With ASV benchmarks, users can pick a few representative commits from `main` to build a performance timeline. For example, use `--steps 3` to sample three commits evenly spaced along the branch. ASV will use the branch(es) listed in `asv.conf.json` (currently `main`) to generate the summary graphs.

```bash
cd benchmarks
git checkout main
uv run asv run main --steps 3
uv run asv publish
uv run asv preview --browser
```

> [!WARNING]
> A full run executes every benchmark on Python 3.11 and 3.14 and can take some time! Use the `--quick` flag to run a subset of benchmarks.

### Compare two refs

`asv continuous` benchmarks two git refs side by side and prints a summary of benchmarks that are faster, slower, or unchanged. `main` is the baseline; the second argument is the ref to test.

With your PR branch checked out:

```bash
uv run asv continuous main HEAD
```

Or compare `main` against a named branch without checking it out:

```bash
uv run asv continuous main <new-branch>
```

In the example above, ASV benchmarks `<new-branch>` and compares it against main as the
baseline.

In both cases, `HEAD` when your PR branch is checked out or name the branch explicitly, the
`main` branch is the baseline and the second reference is what gets benchmarked against it.

For a quicker but less trustworthy check:

```bash
uv run asv continuous main HEAD --quick
```

Useful flags:

```bash
uv run asv run --quick -e                  # fast check; results not saved
uv run asv run --quick --show-stderr       # fast smoke check; results not saved
uv run asv run --show-stderr -b Small      # debug failures (regex on benchmark name)
```

## Clean previous outputs

Delete generated ASV artifacts when environments go stale, results look inconsistent, or you want a fresh local run:

```bash
cd benchmarks
rm -rf env results html virtualenv .asv
```

This removes benchmark environments (`env/`), saved timings (`results/`), and the published site (`html/`). ASV recreates them on the next run. The `virtualenv/` directory is a leftover from older configs.

If you removed `.asv/`, you can re-register your machine before benchmarking again:

```bash
uv run asv machine --yes
```

Note that all of these paths are gitignored. See [Output layout](#output-layout) for what each directory contains.

## Output layout

- `env/`: ASV uv virtualenvs.
- `results/`: Raw JSON timing data.
- `html/`: Published site (`asv publish`).
- `virtualenv/`: Leftover from older ASV configs.
- `.asv/`: Local machine metadata.

## Benchmark modules

- `MembraneCurvatureSmallBenchmark`: `MembraneCurvature.run()` on a nine-atom PO4 system (Fourier and binning).
- `MembraneCurvatureBenchmark`: `MembraneCurvature.run()` on ~900 PO4 atoms at grid sizes 25, 50, and 100.
- `MembraneCurvatureFourierModesBenchmark`: Fourier mode truncation `(2, 2)`, `(3, 3)`, and `(5, 5)`.
- `MembraneCurvatureTrajectoryBenchmark`: Multi-frame Fourier and binning runs, including peak memory.

The file `membranecurvature.py` contains the benchmarks for the `MembraneCurvature.run()` method for different input sizes and for the two surface derivation methods. Peak memory benchmarks are also included to track RAM usage during multi-frame trajectory runs, where per-frame surface and curvature arrays grow with trajectory length and grid resolution.

## Pull requests

### Before opening or updating a PR

1. **Check for correctness (fast, required for benchmark changes)**
   Run the benchmarks with the `--quick` flag to check for correctness:

   ```bash
   uv run asv run --quick --show-stderr
   ```

   With `--quick`, results are not saved. The run only checks that every benchmark completes.

2. **Optional: full suite smoke on your branch**

   If you changed benchmark code or expect a performance impact, benchmark the current commit:

   ```bash
   uv run asv run HEAD --steps 1 --show-stderr
   ```

   This runs the full suite on Python 3.11 and 3.14 and saves one commit's timings. **It is not required on every PR, but it is recommended to run it if you changed benchmark code or expect a performance impact.**

3. **Optional: performance comparison vs `main`**

   When the PR may change speed, compare your branch against `main` (see [Compare two refs](#compare-two-refs)):

   ```bash
   uv run asv continuous main HEAD
   ```

   ASV compares `HEAD` (your branch) against `main` and prints a summary table of benchmarks that are faster, slower, or unchanged. This runs locally and does **not** modify the public `main` history. This is a quick comparison to be included in the PR discussion (optional).


[uv]:https://docs.astral.sh/uv/
[uv plugin]:https://asv.readthedocs.io/en/latest/autoapi/asv/plugins/uv/index.html
[Airspeed Velocity]:https://asv.readthedocs.io/
[asv.conf.json]: asv.conf.json
[pyproject.toml]:../pyproject.toml
