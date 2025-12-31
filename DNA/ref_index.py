from Utils.utils import run_quiet, _resolve_outdir
from pathlib import Path
import shutil

class RefIndex: 
    def __init__(self, args):
        self.args = args
        self.outdir = _resolve_outdir(base_outdir=args.outdir)
        self.threads = args.threads
        self.ref_fasta = Path(args.ref_fasta).resolve()
        

    def run(self):
        ref_index_outdir = _resolve_outdir(
            base_outdir=self.outdir,
            subdir="01_DNAseq_alignment",
            ensure_dir=True,
        )

        # Write indices in the pipeline output directory to avoid permission issues on read-only refs
        dest_ref = ref_index_outdir / self.ref_fasta.name
        if dest_ref.resolve() != self.ref_fasta.resolve():
            if not dest_ref.exists():
                shutil.copy2(self.ref_fasta, dest_ref)
            self.ref_fasta = dest_ref

        # Keep the caller in sync so downstream steps (alignment, DeepVariant, etc.)
        # use the copied and indexed reference path rather than the original.
        self.args.ref_fasta = self.ref_fasta

        cmd1 = f"bwa-mem2 index {self.ref_fasta}"
        cmd2 = f"samtools faidx {self.ref_fasta}"

        run_quiet(cmd1, step="ref_index", logdir=ref_index_outdir)
        run_quiet(cmd2, step="samtools_faidx", logdir=ref_index_outdir)
