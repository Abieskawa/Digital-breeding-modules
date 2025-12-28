# Utils/utils.py
import os
import sys
from typing import Dict
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path  # <-- add this

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
