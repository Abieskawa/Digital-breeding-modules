# Utils/utils.py
import gzip
import os
import numpy as np
import shutil
import subprocess as sbp
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path  # <-- add this
from typing import Dict, Iterable, List, Optional, TYPE_CHECKING, Tuple
from zoneinfo import ZoneInfo

if TYPE_CHECKING:
    import pandas as pd

def time_stamp(msg, logfile=None):
    tzname = os.environ.get("TZ", "UTC")
    try:
        tz = ZoneInfo(tzname)
    except Exception:
        tz = ZoneInfo("UTC")
    stamp = '[Digital breeding pipeline:{}]'.format(
        datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')
    )
    full_msg = f"{stamp} {msg}"
    print(full_msg, flush=True)
    if logfile:
        with open(logfile, "a") as f:
            f.write(full_msg + "\n")
    return full_msg

def read_config(config_file) -> Dict[str, str]:
    """
    Simple key=value loader with support for blank lines, full-line comments,
    and inline comments after '#'.
    """
    cfg: Dict[str, str] = {}
    with open(config_file, "r") as cf:
        for raw in cf:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "#" in line:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg

def parse_int_list(raw: str) -> List[int]:
    if raw is None:
        return []
    s = str(raw).strip()
    if not s:
        return []
    parts = []
    for chunk in s.replace("\t", " ").replace("\n", " ").split(","):
        if chunk.strip():
            parts.extend([p for p in chunk.split(" ") if p.strip()])
    return [int(p) for p in parts]

def _open_text_maybe_gz(path: str):
    """
    Open a text file that may be plain or gzipped, based on magic header.
    """
    with open(path, "rb") as f:
        sig = f.read(2)
    opener = gzip.open if sig == b"\x1f\x8b" else open
    return opener(path, "rt")

def _save_numeric_npz(df_numeric, out_npz: str) -> None:
    os.makedirs(os.path.dirname(out_npz) or ".", exist_ok=True)
    np.savez_compressed(
        out_npz,
        genotype=df_numeric.to_numpy(dtype=np.float32, copy=False),
        samples=np.array(df_numeric.index.astype(str)),
        snps=np.array(df_numeric.columns.astype(str)),
    )

def coerce_bool(raw, default: bool) -> bool:
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    s = str(raw).strip().lower()
    if s == "":
        return default
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "n", "off"}:
        return False
    return default

def extract_vcf_samples(vcf_path: Path) -> List[str]:
    """
    Read only the VCF header to obtain sample IDs (columns 10+).
    Works for plain VCF and bgzip/gz VCF.
    """
    with open(vcf_path, "rb") as f:
        sig = f.read(2)
    opener = gzip.open if sig == b"\x1f\x8b" else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                return cols[9:]
    return []


def _ensure_bgzip_and_index(vcf_path: Path) -> Path:
    """
    Ensure <vcf>.gz exists and is indexed (csi/tbi).
    Returns the gz path.
    """
    gz_path = Path(str(vcf_path) + ".gz")
    if not gz_path.exists():
        time_stamp(f"Creating bgzip VCF: {gz_path}")
        sbp.run(f"bgzip -c {vcf_path} > {gz_path}", shell=True, check=True)
    if not (Path(str(gz_path) + ".csi").exists() or Path(str(gz_path) + ".tbi").exists()):
        time_stamp(f"Indexing bgzip VCF: {gz_path}")
        sbp.run(f"bcftools index {gz_path}", shell=True, check=True)
    return gz_path

def setup_local_env(outdir: Path) -> Tuple[Path, Path]:
    """
    Keep temp/cache writes inside the configured output directory.
    """
    outdir = Path(outdir)
    tmp_dir = outdir / ".tmp"
    mpl_dir = outdir / ".cache" / "matplotlib"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    mpl_dir.mkdir(parents=True, exist_ok=True)
    for key in ("TMPDIR", "TEMP", "TMP"):
        os.environ[key] = str(tmp_dir)
    os.environ["MPLCONFIGDIR"] = str(mpl_dir)
    os.environ.setdefault("MPLBACKEND", "Agg")
    return tmp_dir, mpl_dir

def parallel_process(items: List[tuple], process_func, max_workers: int = 8, raise_on_error: bool = False):
    """
    Perform parallel processing using ProcessPoolExecutor.
    Each item in `items` is a tuple of arguments for `process_func`.
    """
    print(f"Starting parallel processing with {max_workers} workers")
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_func, item) for item in items]
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                print(f"Task failed with exception: {e}")
                if raise_on_error:
                    raise
    print("Parallel processing complete")
    return results

def _resolve_outdir(
    config: Optional[Dict[str, str]] = None,
    key: str = "output_dir",
    base_outdir: Optional[Path] = None,
    default: str = "DGBreeding",
    subdir: Optional[str] = None,
    resolve: bool = False,
    ensure_dir: bool = False,
) -> Path:
    """
    Resolve output directories consistently (supports base, keyed, and subdir paths).
    Optionally create the directory when ensure_dir=True.
    """
    if base_outdir is None:
        if config is None:
            base_outdir = Path(default)
        else:
            if key == "output_dir":
                raw = config.get(key, None)
                raw = str(raw).strip() if raw is not None else ""
                base_outdir = Path(raw or default)
            else:
                raw_base = config.get("output_dir", None)
                raw_base = str(raw_base).strip() if raw_base is not None else ""
                base_outdir = Path(raw_base or default)
    else:
        base_outdir = Path(base_outdir)

    base_outdir = base_outdir.expanduser()

    if config is not None and key != "output_dir":
        raw = config.get(key, None)
        raw = str(raw).strip() if raw is not None else ""
        if raw:
            override = Path(raw).expanduser()
            base_outdir = override if override.is_absolute() else base_outdir / override

    if subdir:
        base_outdir = base_outdir / subdir

    outdir = base_outdir.resolve() if resolve else base_outdir
    if ensure_dir:
        outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def run_quiet(cmd, workdir=None, step='cmd', logdir=None):
    """Run a shell command quietly, optionally logging stdout/stderr to a file."""
    import subprocess as sbp

    cmd = " ".join(str(cmd).split())
    log_handle = None
    stdout = sbp.DEVNULL
    stderr = sbp.DEVNULL

    if logdir:
        log_path_parent = Path(logdir)
        log_path_parent.mkdir(parents=True, exist_ok=True)
        log_path = log_path_parent / f"{step}.log"
        # open in binary to safely capture arbitrary subprocess output
        log_handle = open(log_path, 'wb')
        stdout = log_handle
        stderr = sbp.STDOUT

    try:
        sbp.run(cmd,
                shell=True,
                check=True,
                stdout=stdout,
                stderr=stderr,
                cwd=str(workdir) if workdir else None)
    except sbp.CalledProcessError:
        if logdir:
            call_log(logdir, step, cmd)
        sys.exit(1)
    finally:
        if log_handle:
            log_handle.close()

def _append_log_line(log_dir: Path, step: str, msg: str) -> Path:
    """
    Append a human-readable line to <log_dir>/<step>.log (creating dirs as needed).
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{step}.log"
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(msg.rstrip() + "\n")
    return log_path

def call_log(out_dir, name, cmd):
    # always give time_stamp a message
    cmd_clean = " ".join(str(cmd).split())
    time_stamp(f"!!ERROR!! Command failed ({name}): {cmd_clean}")

    # try the current step layout first: <out_dir>/<name>.log
    out_dir = Path(out_dir)
    candidates = [
        out_dir / f"{name}.log",
        out_dir / "log" / f"{name}.log",  # legacy fallback
    ]
    log_path = next((p for p in candidates if p.exists()), None)

    print('please check below:\n')
    if not log_path:
        print(f"(no log found at: {candidates[0]} or {candidates[1] if len(candidates)>1 else ''})")
        return

    # tail the log (avoid huge dumps)
    try:
        with open(log_path, "r", errors="replace") as log:
            for line in log:
                print(line, end='')
    except Exception as e:
        print(f"(warn) could not read log {log_path}: {e}")


# --------------------------- FASTQ discovery ---------------------------

# Cleaned outputs from fastp (priority)
R1_PATTERNS_CLEANED = [
    "*_R1.cleaned.gz",
    "*_1.cleaned.gz",
    "*_R1.cleaned.fastq.gz",
    "*_1.cleaned.fastq.gz",
]

# Common raw patterns (fallback)
R1_PATTERNS_RAW = [
    "*_1.fastq", "*_1.fq", "*_1.fastq.gz", "*_1.fq.gz",
    "*_R1.fastq", "*_R1.fq", "*_R1.fastq.gz", "*_R1.fq.gz",
]


def _discover_pairs(
    wd: Path,
    r1_patterns: List[str],
    *,
    require_r2: bool,
    sort_paths: bool = True,
) -> List[Tuple[str, str, Path, Optional[Path]]]:
    wd = Path(wd)
    seen = set()
    pairs: List[Tuple[str, str, Path, Optional[Path]]] = []
    for pat in r1_patterns:
        it = sorted(wd.glob(pat)) if sort_paths else wd.glob(pat)
        for r1p in it:
            fname = r1p.name
            if "_R1." in fname:
                base, ext = fname.split("_R1.", 1)
                r2n = f"{base}_R2.{ext}"
                style = "R"
            elif "_1." in fname:
                base, ext = fname.split("_1.", 1)
                r2n = f"{base}_2.{ext}"
                style = "1"
            else:
                continue
            r2p = wd / r2n
            if require_r2 and not r2p.exists():
                print(f"[info] Skip {fname}: missing pair {r2n}", file=sys.stderr)
                continue
            key = (base, r1p.name, r2p.name if r2p.exists() else "")
            if key in seen:
                continue
            seen.add(key)
            pairs.append((base, style, r1p, r2p if r2p.exists() else None))
    return pairs


def _discover_pairs_in_dir(wd: Path) -> List[Dict]:
    """
    Return list of {sample, r1, r2|None}, preferring fastp-cleaned names,
    otherwise falling back to raw naming patterns.
    'sample' is the basename before _R1/_1 suffix + extension.
    """
    wd = Path(wd)
    pairs: List[Tuple[str, str, Path, Optional[Path]]] = []

    for patterns in (R1_PATTERNS_CLEANED, R1_PATTERNS_RAW):
        pairs = _discover_pairs(wd, patterns, require_r2=False, sort_paths=True)
        if pairs:
            break

    if pairs:
        return [dict(sample=base, r1=r1p, r2=r2p) for base, _, r1p, r2p in pairs]

    all_fqs = (
        sorted(wd.glob("*.fastq"))
        + sorted(wd.glob("*.fq"))
        + sorted(wd.glob("*.fastq.gz"))
        + sorted(wd.glob("*.fq.gz"))
    )
    return [dict(sample=p.stem, r1=p, r2=None) for p in all_fqs]


def discover_pairs_by_patterns(wd: Path, r1_patterns: List[str]) -> List[Tuple[str, str, str, str]]:
    """
    Find paired-end files in wd, given R1 glob patterns.

    Returns list of (basename, style, r1_name, r2_name)
      style = "R" if filenames use _R1/_R2, else "1" for _1/_2
    """
    pairs = _discover_pairs(wd, r1_patterns, require_r2=True, sort_paths=False)
    results: List[Tuple[str, str, str, str]] = []
    for base, style, r1p, r2p in pairs:
        if r2p is None:
            continue
        results.append((base, style, r1p.name, r2p.name))
    return results


def _kreport_total_reads(kreport: Path) -> Optional[int]:
    """
    Return the total number of sequences Kraken2 saw (root reads_clade from kreport).
    """
    try:
        with open(kreport) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    total = int(parts[1])
                except ValueError:
                    continue
                return total
    except Exception as e:
        time_stamp(f"[contam] failed to read kreport {kreport}: {e}")
    return None


# --------------------------- Plot helpers ---------------------------

def _annotate_box_fliers(ax, x_values: Iterable, y_values: Iterable, labels: Iterable, rotation: int = 45) -> None:
    """
    Label boxplot fliers on the given axes.
    x_values / y_values / labels must be aligned and of equal length.
    """
    try:
        import numpy as np
        from matplotlib.cbook import boxplot_stats
    except Exception:
        return

    try:
        from adjustText import adjust_text
    except Exception:
        adjust_text = None

    x_arr = np.asarray(list(x_values))
    y_arr = np.asarray(list(y_values), dtype=float)
    labels_arr = np.asarray(list(labels))

    if x_arr.size == 0 or y_arr.size == 0 or labels_arr.size == 0:
        return

    texts = []
    for x_val in np.unique(x_arr):
        mask = x_arr == x_val
        if not np.any(mask):
            continue
        x_vals = x_arr[mask]
        y_vals = y_arr[mask]
        labs = labels_arr[mask]
        bstats = boxplot_stats(y_vals, labels=None)[0]
        fliers = bstats.get("fliers", [])
        if not fliers:
            continue
        fliers_set = set(fliers)
        for xv, yv, lab in zip(x_vals, y_vals, labs):
            if yv in fliers_set:
                texts.append(ax.text(xv, yv, lab, rotation=rotation, ha="right", va="center", fontsize=7))

    if texts and adjust_text:
        try:
            adjust_text(
                texts,
                arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5),
                force_points=(0.1, 0.2),
                force_text=(0.05, 0.2),
                expand_points=(1.05, 1.1),
                lim=150,
            )
        except Exception:
            pass


def _count_reads_in_fastq(path: Path) -> int:
    """Number of reads = lines/4 (handles .gz, uses pigz for speed)."""
    if str(path).endswith(".gz"):
        pigz = shutil.which("pigz")
        if not pigz:
            time_stamp("[contam] pigz not found on PATH; install pigz in the environment.")
            return 0
        cmd = f"{pigz} -cd {path} | wc -l"
    else:
        cmd = f"wc -l {path} | awk '{{print $1}}'"
    try:
        out = sbp.check_output(cmd, shell=True, text=True).strip()
        nlines = int(out.split()[0])
        return nlines // 4
    except Exception as e:
        time_stamp(f"[contam] failed to count reads in {path}: {e}")
        return 0
