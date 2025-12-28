#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path


class MappingEvalBase:
    def _setup_local_env(self) -> None:
        """
        Keep temp/cache writes inside the configured output directory.
        """
        tmp_dir = self.outdir / ".tmp"
        mpl_dir = self.outdir / ".cache" / "matplotlib"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        mpl_dir.mkdir(parents=True, exist_ok=True)
        for key in ("TMPDIR", "TEMP", "TMP"):
            os.environ[key] = str(tmp_dir)
        os.environ["MPLCONFIGDIR"] = str(mpl_dir)
        os.environ.setdefault("MPLBACKEND", "Agg")
        self.tmp_dir = tmp_dir
        self.mpl_config_dir = mpl_dir

    def __init__(self, args):
        self.threads = int(getattr(args, "threads", 1))
        self.outdir = Path(getattr(args, "outdir", ".")).resolve()
        self.input_dir = Path(getattr(args, "input_dir", ".")).resolve()

        self._setup_local_env()

        self.fastp_outdir = Path(getattr(args, "fastp_outdir", self.input_dir)).resolve()
        self.align_dir = Path(getattr(args, "align_outdir", self.outdir)).resolve()

        self.gff = Path(getattr(args, "gff", "")).expanduser().resolve() if getattr(args, "gff", None) else None

        raw_tag = getattr(args, "tag_outliers", True)
        if isinstance(raw_tag, str):
            self.tag_outliers = raw_tag.strip().lower() not in {"0", "false", "no", "off", ""}
        else:
            self.tag_outliers = bool(raw_tag)

        default_kraken_db = "/kraken_db"
        db_val = getattr(args, "kraken_db", default_kraken_db)
        if isinstance(db_val, str):
            db_val = (db_val or default_kraken_db).strip()
        else:
            db_val = str(db_val).strip() if db_val else default_kraken_db
        self.kraken_db = db_val

        self.report_dir = self.outdir / "evaluation"
        self.report_dir.mkdir(parents=True, exist_ok=True)
        (self.report_dir / "logs").mkdir(parents=True, exist_ok=True)
