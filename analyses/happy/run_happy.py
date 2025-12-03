#!/usr/bin/env python3

import subprocess
import os
import hydra
from omegaconf import DictConfig


def run_happy(truth_vcf, truth_bed, query_vcf,
              stratifications, output_prefix,
              reference_fasta, sdf, threads=6):
    """Run hap.py with vcfeval engine."""
    hap_env = 'hap.py'
    cmd = [
        "mamba", "run", "-n", hap_env, "hap.py",
        truth_vcf,
        query_vcf,
        "-f", truth_bed,
        "-r", reference_fasta,
        "-o", output_prefix,
        "-t", "ga4gh",
        "-V",
        "--pass-only",
        "--fixchr",
        "--stratification", stratifications,
        "--engine", "vcfeval",
        "--engine-vcfeval-template", sdf,
        "--threads", str(threads),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


@hydra.main(config_path="conf", config_name="config", version_base="1.3")
def main(cfg: DictConfig):
    reference = cfg.reference
    threads = cfg.threads
    stratifications = cfg.stratifications
    sdf = cfg.sdf

    for sample, paths in cfg.samples.items():
        truth_vcf = paths.truth_vcf
        truth_bed = paths.truth_bed
        query_vcf = paths.query_vcf
        output_prefix = os.path.join(cfg.output_dir, sample)

        os.makedirs(cfg.output_dir, exist_ok=True)
        run_happy(truth_vcf, truth_bed, query_vcf,
                  stratifications, output_prefix,
                  reference, sdf, threads)


if __name__ == "__main__":
    main()
