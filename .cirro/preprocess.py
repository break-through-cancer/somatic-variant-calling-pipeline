"""
preprocess.py -- Cirro preprocess hook for the somatic variant calling pipeline.

Runs on Cirro before Nextflow starts. Reads the input dataset's file manifest (ds.files:
one row per file, with `sample` and `file` columns) and injects the parameters main.nf needs:

    runs                 [{sample_name, tumor_bam, tumor_bam_index}, ...]   one per tumor
    normal_sample_name   name of the shared normal
    normal_bam           path of the shared normal BAM/CRAM
    normal_bam_index     its index
    callers              comma-separated callers picked on the launch form

Input: preprocessed (duplicate-marked, recalibrated) BAM or CRAM files, each with its index.
Shared-normal model: every tumor is paired with the SAME single normal. Any sample whose name
contains one of NORMAL_NAME_MARKERS is the normal; everything else is a tumor.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

NORMAL_NAME_MARKERS = ["PBMC"]  # adjust if this cohort's normals aren't PBMC-labeled

# Callers that can be ticked on the launch form (use_<caller>), in canonical order.
ALL_CALLERS = ["mutect", "deepsomatic", "strelka"]

# (data extension, index extension), in order of preference: BAM first, then CRAM.
# MuTect1 can only read BAM; the other callers take either.
ALIGNMENT_FORMATS = [(".bam", ".bai"), (".cram", ".crai")]

# Within one format, prefer recalibrated output, then sorted/mapped, then markduplicates.
NAME_PREFERENCE = [".recal.", ".sorted.", ".md."]


def as_bool(value, default=True):
    """Cirro may hand booleans over as real bools or as 'true'/'false' strings."""
    if value is None or value == "":
        return default
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes", "on")
    return bool(value)


def selected_callers(params) -> list:
    return [c for c in ALL_CALLERS if as_bool(params.get(f"use_{c}"))]


def _rank(path):
    name = path.rsplit("/", 1)[-1]
    for i, marker in enumerate(NAME_PREFERENCE):
        if marker in name:
            return i
    return len(NAME_PREFERENCE)


def extract_alignments(ds):
    """Pick one alignment + matching index per sample ({sample: {"bam": ..., "bai": ...}})."""
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    alignments = {}
    for sample, group in df.groupby("sample"):
        files = set(group["file"])
        for data_ext, index_ext in ALIGNMENT_FORMATS:
            # Only accept a data file whose own index sits next to it
            pairs = [f for f in files if f.endswith(data_ext) and f + index_ext in files]
            if pairs:
                best = sorted(pairs, key=lambda f: (_rank(f), f))[0]
                alignments[str(sample)] = {"bam": best, "bai": best + index_ext}
                break

    if not alignments:
        raise ValueError("No BAM/CRAM files with a matching index found in the dataset")
    return alignments


def is_normal(sample_name):
    upper = sample_name.upper()
    return any(marker in upper for marker in NORMAL_NAME_MARKERS)


def main():
    ds = PreprocessDataset.from_running()

    callers = selected_callers(ds.params)
    if not callers:
        raise ValueError("No callers selected; tick at least one of MuTect1, DeepSomatic, Strelka")

    if ds.files.empty:
        raise ValueError(
            "No files found in the input dataset. Check that the dataset's files are "
            "registered to samples in Cirro (More Info > Samplesheet and Files)."
        )

    alignments = extract_alignments(ds)

    normal = None
    tumors = []
    for sample, files in alignments.items():
        if is_normal(sample):
            if normal is not None:
                raise ValueError(
                    f"Multiple normal-labeled samples found ({normal['sample_name']!r} and "
                    f"{sample!r}) -- this pipeline expects exactly one shared normal per "
                    f"dataset. Adjust NORMAL_NAME_MARKERS or the dataset contents."
                )
            normal = {"sample_name": sample, **files}
        else:
            tumors.append({
                "sample_name": sample,
                "tumor_bam": files["bam"],
                "tumor_bam_index": files["bai"],
            })

    if not tumors:
        raise ValueError("No tumor BAM/CRAM found")
    if normal is None:
        raise ValueError(
            "No normal-labeled sample found -- this pipeline requires a matched normal shared "
            "across all tumor samples (see NORMAL_NAME_MARKERS)"
        )

    # MuTect1 reads BAM only: fail now, with the offending files, instead of mid-run.
    if "mutect" in callers:
        cram = [t["tumor_bam"] for t in tumors if t["tumor_bam"].endswith(".cram")]
        if normal["bam"].endswith(".cram"):
            cram.append(normal["bam"])
        if cram:
            raise ValueError(
                "MuTect1 cannot read CRAM, but these inputs are CRAM: "
                + ", ".join(cram[:3]) + (" ..." if len(cram) > 3 else "")
                + ". Untick MuTect1 or provide BAMs."
            )

    ds.add_param("runs", tumors)
    ds.add_param("normal_sample_name", normal["sample_name"])
    ds.add_param("normal_bam", normal["bam"])
    ds.add_param("normal_bam_index", normal["bai"])
    ds.add_param("callers", ",".join(callers))

    print(f"[INFO] Callers: {', '.join(callers)}")
    print(f"[INFO] Normal: {normal['sample_name']}")
    print(f"[INFO] Tumors ({len(tumors)}): {', '.join(t['sample_name'] for t in tumors)}")


if __name__ == "__main__":
    main()
