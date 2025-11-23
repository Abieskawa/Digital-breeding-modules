import sys
import subprocess as sbp # enable activate and terminate program
from Utils.utils import clean_cmd, call_log
from pathlib import Path
import shutil

class RefIndex: 
    def __init__(self, args):
        self.args = args
        self.outdir = Path(args.outdir)
        self.threads = args.threads
        self.ref_fasta = Path(args.ref_fasta).resolve()
        

    def run(self):
        ref_index_outdir = Path(self.outdir / "01_DNAseq_alignment")

        if not ref_index_outdir.exists():
            ref_index_outdir.mkdir(parents=True, exist_ok=True)

        ref_index_log = ref_index_outdir / "ref_index.log"
        faidx_log = ref_index_outdir / "samtools_faidx.log"

        # Write indices in the pipeline output directory to avoid permission issues on read-only refs
        dest_ref = ref_index_outdir / self.ref_fasta.name
        if dest_ref.resolve() != self.ref_fasta.resolve():
            if not dest_ref.exists():
                shutil.copy2(self.ref_fasta, dest_ref)
            self.ref_fasta = dest_ref

        # Keep the caller in sync so downstream steps (alignment, DeepVariant, etc.)
        # use the copied and indexed reference path rather than the original.
        self.args.ref_fasta = self.ref_fasta

        cmd1 = (
            f"bwa-mem2 index {self.ref_fasta} >> {ref_index_log} 2>&1"
        )

        cmd2 = (
            f"samtools faidx {self.ref_fasta} >> {faidx_log} 2>&1"
        )

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(ref_index_outdir, 'ref_index', cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(ref_index_outdir, 'samtools_faidx', cmd2)
            sys.exit(1)
