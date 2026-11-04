# somatic-variant-calling-pipeline

Somatic variant calling on preprocessed tumor/normal BAM or CRAM files, with an optional
consensus across the callers. One Nextflow DSL2 workflow, packaged for Cirro.

It replaces running the separate `mutect-1`, `deepsomatic`, sarek-Strelka and
`consensus-calling` pipelines and stitching their outputs together by hand: everything is
discovered from **one input dataset**, so no master samplesheet is needed.

```
BAM/CRAM  ->  MuTect1 ----------------\
          ->  DeepSomatic -------------+->  consensus (optional)
          ->  Manta -> Strelka2 -------/
```

## Input

One Cirro dataset of **preprocessed** (duplicate-marked, recalibrated) BAM or CRAM files, each
with its index. Alignment is not done here.

- Any sample whose name contains `PBMC` is the **normal**; everything else is a tumor
  (`NORMAL_NAME_MARKERS` in `.cirro/preprocess.py`).
- There must be exactly **one shared normal**; every tumor is paired with it.
- BAM is preferred over CRAM when a sample has both. **MuTect1 needs BAM**; the other callers take
  either, and the run fails at preprocess if MuTect1 is ticked and the inputs are CRAM.

`.cirro/preprocess.py` turns `ds.files` into the params `main.nf` reads (`runs`,
`normal_sample_name`, `normal_bam`, `normal_bam_index`, `callers`).

## Callers

| Caller | Reports | Notes |
|---|---|---|
| MuTect1 | SNVs | scatter/gather + optional ContEst; code from the standalone `mutect-1` pipeline |
| DeepSomatic | SNVs, indels | NVIDIA Parabricks, needs a GPU; code from the standalone `deepsomatic` pipeline |
| Manta -> Strelka2 | SNVs, indels (+ SVs from Manta) | Manta runs first and its candidate small indels seed Strelka (`--indelCandidates`), as in nf-core/sarek |

Pick any combination on the launch form. One optional set of calling intervals (BED or interval list, `target_list`) applies to MuTect1, Manta and
Strelka (converted/compressed/indexed automatically for Manta/Strelka); `wes` switches Manta/Strelka to targeted mode.

## Consensus

On when **Generate consensus** is ticked and **two or more callers** are selected. It runs on top
of each caller's own results, which are always published.

How callers are combined depends on how many take part (k):

| k | rule |
|---|---|
| 3 | a variant must be called by at least 2 callers (n-1) |
| 2 | both callers must call it (intersection) |
| 1 | the caller's calls pass through (only happens for indels, e.g. MuTect1 + Strelka) |

SNVs use every selected caller. Indels use only the callers that report indels (not MuTect1).
Only `PASS` records count by default (`pass_only`), so DeepSomatic's `RefCall` records don't
count as support. The final VCF is the SNV consensus plus the indel consensus.

## Outputs

Under `<outdir>/<sample>/`:

- `<sample>.mutect1.pass.vcf`, `<sample>.mutectv1.final.vcf`, call stats (MuTect1)
- `<sample>.deepsomatic.vcf.gz` (DeepSomatic)
- `<sample>_vs_<normal>.strelka.somatic_snvs.vcf.gz` / `.somatic_indels.vcf.gz` (Strelka)
- `manta/<sample>.manta.{somatic_sv,diploid_sv,candidate_sv,candidate_small_indels}.vcf.gz`
- `<sample>.final_consensus.vcf.gz`, `<sample>.consensus_log.txt`, `snv_*` / `indel_*` consensus files

The caller names are in the filenames, the same as the standalone pipelines, so existing
downstream tools keep working.

## Layout

```
main.nf                         workflow: validation, caller selection, consensus
nextflow.config                 param defaults
modules/local/                  mutect1.nf, deepsomatic.nf, manta.nf, strelka.nf, consensus.nf
subworkflows/local/             mutect1_calling, manta_strelka, consensus_calling
.cirro/                         preprocess.py, process-form.json, process-input.json, process-compute.config
```

## Registering on Cirro

The process needs `usesSampleSheet: true`, and its parent process types must include whatever
produces your BAM/CRAM dataset (e.g. `process-nf-core-sarek-align-3-2`, `jenc-recalibrate`,
`ingest_aligned_bam`) so those datasets can be selected as input. Those producers need file
mapping rules that register the BAM/CRAM and index files to samples; otherwise `ds.files` is empty.

## Not included from the standalone pipelines

- MuTect1's force-calling-interval helpers (`make_force_calling_intervals`).
- The duplicate-marking check (inputs are assumed to be preprocessed).
- Alignment, BQSR and anything else upstream of variant calling.
