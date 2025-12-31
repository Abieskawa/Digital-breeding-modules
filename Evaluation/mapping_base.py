#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from Utils.utils import _resolve_outdir, setup_local_env


class MappingEvalBase:
    def _setup_local_env(self) -> None:
        self.tmp_dir, self.mpl_config_dir = setup_local_env(self.outdir)

    def __init__(self, args):
        self.threads = int(getattr(args, "threads", 1))
        self.outdir = _resolve_outdir(base_outdir=getattr(args, "outdir", "."), resolve=True)
        self.input_dir = Path(getattr(args, "input_dir", ".")).resolve()

        self._setup_local_env()

        self.fastp_outdir = _resolve_outdir(base_outdir=getattr(args, "fastp_outdir", self.input_dir), resolve=True)
        self.align_dir = _resolve_outdir(base_outdir=getattr(args, "align_outdir", self.outdir), resolve=True)

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

        self.report_dir = _resolve_outdir(
            base_outdir=self.outdir,
            subdir="evaluation",
            resolve=True,
            ensure_dir=True,
        )
        (self.report_dir / "logs").mkdir(parents=True, exist_ok=True)
