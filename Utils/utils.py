# Utils/utils.py
import os
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

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def call_log(out_dir, name, cmd):
    # always give time_stamp a message
    time_stamp(f"!!ERROR!! Command failed ({name}): {clean_cmd(cmd)}")

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
