#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compatibility wrapper for legacy MappingEvaluator usage.
"""
from Evaluation.contamination import ContaminationEvaluator
from Evaluation.fastq_qc import FastQCEvaluator
from Evaluation.mapping_yield import MappingYieldEvaluator

__all__ = ["MappingEvaluator"]


class MappingEvaluator:
    """
    Steps:
      0: FastQC + MultiQC
      1: Kraken2 on all reads; aggregate species Top-5 across samples (boxplot)
      2: Mapping yield ratio = unique (primary, non-supplementary) BAM reads / fastp reads; one boxplot
    """

    def __init__(self, args):
        self._qc = FastQCEvaluator(args)
        self._contam = ContaminationEvaluator(args)
        self._align = MappingYieldEvaluator(args)

    def evaluate_qc(self, fastq_stage: str = "processed") -> None:
        self._qc.evaluate_qc(fastq_stage=fastq_stage)

    def evaluate_contamination(self):
        return self._contam.evaluate_contamination()

    def evaluate_alignment(self):
        return self._align.evaluate_alignment()
