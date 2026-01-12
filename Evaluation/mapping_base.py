#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from Utils.utils import _resolve_outdir, setup_eval_env

DEFAULT_KRAKEN_DB = "/kraken2_db"


class MappingEvalBase:
    def __init__(self, args):
        self.threads = int(getattr(args, "threads", 1))
        self.outdir = _resolve_outdir(base_outdir=getattr(args, "outdir", "."), resolve=True)
        self.input_dir = Path(getattr(args, "input_dir", ".")).resolve()

        self.tmp_dir, self.mpl_config_dir, self.report_dir = setup_eval_env(self.outdir)

        self.fastp_outdir = _resolve_outdir(base_outdir=getattr(args, "fastp_outdir", self.input_dir), resolve=True)
        self.align_dir = _resolve_outdir(base_outdir=getattr(args, "align_outdir", self.outdir), resolve=True)

        self.gff = Path(getattr(args, "gff", "")).expanduser().resolve() if getattr(args, "gff", None) else None

        raw_tag = getattr(args, "tag_outliers", True)
        if isinstance(raw_tag, str):
            self.tag_outliers = raw_tag.strip().lower() not in {"0", "false", "no", "off", ""}
        else:
            self.tag_outliers = bool(raw_tag)

        self.kraken_db = DEFAULT_KRAKEN_DB
        # report_dir already created by setup_eval_env
