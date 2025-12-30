#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess as sbp
import sys
from pathlib import Path
from typing import List

from Utils.utils import _discover_pairs_in_dir, call_log, time_stamp
from Evaluation.mapping_base import MappingEvalBase


class FastQCEvaluator(MappingEvalBase):
    """
    Steps:
      0/1: FastQC + MultiQC
    """

    def evaluate_qc(self, fastq_stage: str = "processed") -> None:
        """
        Run FastQC + MultiQC on FASTQs from a specified stage.
        fastq_stage:
          - 'raw': only consider the original input directory
          - 'processed': only consider fastp-cleaned outputs
        """
        if fastq_stage is None:
            time_stamp("[qc] fastq_stage is required; expected one of: raw, processed")
            sys.exit(1)
        stage = str(fastq_stage).strip().lower()
        if not stage:
            time_stamp("[qc] fastq_stage is blank; expected one of: raw, processed")
            sys.exit(1)
        choices = {"raw", "processed"}
        if stage not in choices:
            time_stamp(
                f"[qc] invalid fastq_stage '{stage}'; expected one of: {', '.join(sorted(choices))}"
            )
            sys.exit(1)

        stage_dirs = []
        if stage == "processed":
            stage_dirs.append(("processed", self.fastp_outdir))
        if stage == "raw":
            stage_dirs.append(("raw", self.input_dir))

        pairs = []
        for label, wd in stage_dirs:
            pairs = _discover_pairs_in_dir(wd)
            if pairs:
                break

        fastqs: List[Path] = []
        for d in pairs:
            fastqs.append(d["r1"])
            if d["r2"]:
                fastqs.append(d["r2"])
        if not fastqs:
            time_stamp(f"[qc] no FASTQs found for stage '{stage}'; skip")
            return

        stage_dir = stage
        fq_out = self.report_dir / "fastqc" / stage_dir
        mq_out = self.report_dir / "multiqc" / stage_dir
        fq_out.mkdir(parents=True, exist_ok=True)
        mq_out.mkdir(parents=True, exist_ok=True)

        files_str = " ".join(str(p) for p in fastqs)
        cmd1 = f"fastqc -t {self.threads} -o {fq_out} {files_str}"
        time_stamp(f"[qc:{stage_dir}] {cmd1}")
        try:
            sbp.run(cmd1, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(self.report_dir / "logs", f"fastqc_{stage_dir}", cmd1)

        cmd2 = f"multiqc -o {mq_out} {fq_out}"
        time_stamp(f"[qc:{stage_dir}] {cmd2}")
        try:
            sbp.run(cmd2, shell=True, check=True)
        except sbp.CalledProcessError:
            call_log(self.report_dir / "logs", f"multiqc_{stage_dir}", cmd2)
