"""Genomic feature extraction for the AMR_ciprofloxacin isolates (command-line version).

Same method and settings as 02_feature_extraction_colab.ipynb, packaged to run unattended on
a Linux server such as AWS EC2. For each isolate a subset of reads is streamed from the
European Nucleotide Archive (ENA) and three kinds of features are extracted:

  - amino-acid substitutions in gyrA, gyrB, parC, parE, marR, acrR and soxR
    (BWA-MEM against E. coli K-12 MG1655, K-12 numbering)
  - acquired quinolone resistance genes (KMA against the ResFinder database)
  - a Mash sketch for genome distances / lineage clustering

Results are saved per isolate, so the script can be stopped and restarted at any time.
Outputs (in --out): features.csv, qc.csv, mash_distances.csv

Usage:
    bash aws/setup.sh
    nohup ~/venv/bin/python aws/extract_features.py --workers 16 > run.log 2>&1 &
"""
import argparse
import io
import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import requests
from Bio import SeqIO
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Settings (identical to notebook 02)
# ---------------------------------------------------------------------------
TARGET_COVERAGE = 40
MAX_READ_PAIRS = 1_000_000
GENOME_SIZE = 5.0e6

MIN_DEPTH = 5          # minimum reads covering a position
MIN_BASE_FRAC = 0.8    # minimum fraction of reads supporting the majority base
MIN_BASE_QUAL = 20     # ignore bases below this Phred quality
MAX_MISMATCH = 0.10    # discard reads with more than 10% mismatches to the reference gene

MIN_IDENTITY = 90.0    # acquired genes (ResFinder defaults)
MIN_COVERAGE = 60.0

RANDOM_SEED = 42
TARGET_GENES = ["gyrA", "gyrB", "parC", "parE", "marR", "acrR", "soxR"]
FLANK = 150

DATA_URL = "https://raw.githubusercontent.com/Aleroebijuwa/Explainable-ml-ciprofloxacin-resistance/main/Data"
GENBANK_URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               "?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text")
RESFINDER = "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--out", default=str(Path.home() / "amr_cipro"), help="results folder")
parser.add_argument("--workers", type=int, default=os.cpu_count(), help="isolates processed in parallel")
parser.add_argument("--n-isolates", type=int, default=None, help="process only a random subset (default: all)")
parser.add_argument("--no-retry-errors", action="store_true", help="skip isolates that failed previously")
parser.add_argument("--combine-only", action="store_true", help="only build the CSV files from finished isolates")
args = parser.parse_args()

OUT = Path(args.out)
RESULTS_DIR, SKETCH_DIR = OUT / "per_sample", OUT / "sketches"
WORK, REF = OUT / "work", OUT / "reference"
for d in [RESULTS_DIR, SKETCH_DIR, WORK, REF]:
    d.mkdir(parents=True, exist_ok=True)

N_WORKERS = max(1, args.workers)
TOOL_THREADS = 1  # one thread per tool; parallelism comes from running isolates side by side


def log(message):
    print(time.strftime("%H:%M:%S"), message, flush=True)


def sh(cmd):
    # Run a bash command and raise an error with its output if it fails
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\n{result.stderr[-2000:]}")
    return result.stdout


# ---------------------------------------------------------------------------
# Isolates
# ---------------------------------------------------------------------------
metadata = pd.read_csv(f"{DATA_URL}/supplementary_meta_table.csv")
sequencing_metadata = pd.read_csv(f"{DATA_URL}/meta.tsv", sep="\t")
isolates = metadata.merge(sequencing_metadata[["sample_alias", "fastq_ftp"]], on="sample_alias", how="left")


def paired_fastq_urls(fastq_ftp):
    files = str(fastq_ftp).split(";")
    r1 = [f for f in files if f.endswith("_1.fastq.gz")]
    r2 = [f for f in files if f.endswith("_2.fastq.gz")]
    if not r1 or not r2:
        return None, None
    return "https://" + r1[0], "https://" + r2[0]


isolates[["url_r1", "url_r2"]] = isolates["fastq_ftp"].apply(lambda s: pd.Series(paired_fastq_urls(s)))
isolates["run_accession"] = isolates["url_r1"].str.extract(r"/([EDS]RR\d+)_1")[0]
isolates["sample_id"] = isolates["sample_alias"].str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)
isolates = isolates.sample(frac=1, random_state=RANDOM_SEED).reset_index(drop=True)
log(f"Isolates: {len(isolates)}, with paired FASTQ links: {isolates['url_r1'].notna().sum()}")

# ---------------------------------------------------------------------------
# Reference sequences
# ---------------------------------------------------------------------------
mg1655 = SeqIO.read(io.StringIO(requests.get(GENBANK_URL, timeout=300).text), "genbank")
chrom_targets = {}
with open(REF / "chromosomal_targets.fa", "w") as fa:
    for feature in mg1655.features:
        gene = feature.qualifiers.get("gene", [None])[0]
        if feature.type != "CDS" or gene not in TARGET_GENES:
            continue
        start, end = int(feature.location.start), int(feature.location.end)
        region = mg1655.seq[start - FLANK:end + FLANK]
        if feature.location.strand == -1:
            region = region.reverse_complement()
        cds = str(region[FLANK:FLANK + (end - start)])
        assert str(Seq(cds).translate(table=11, cds=True)) == feature.qualifiers["translation"][0]
        chrom_targets[gene] = {"cds": cds, "contig_length": len(region)}
        fa.write(f">{gene}\n{region}\n")
assert set(chrom_targets) == set(TARGET_GENES)


def read_fasta_text(text):
    records, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = []
        elif name:
            records[name].append(line.strip())
    return {k: "".join(v) for k, v in records.items()}


quinolone_db = read_fasta_text(requests.get(f"{RESFINDER}/quinolone.fsa", timeout=300).text)
aminoglycoside_db = read_fasta_text(requests.get(f"{RESFINDER}/aminoglycoside.fsa", timeout=300).text)
aac6_ib = {k: v for k, v in aminoglycoside_db.items() if re.match(r"aac\(6'\)-?Ib(-cr)?_\d", k)}
with open(REF / "acquired_genes.fa", "w") as fa:
    for name, seq in {**quinolone_db, **aac6_ib}.items():
        fa.write(f">{name}\n{seq}\n")
sh(f"cd '{REF}' && bwa index chromosomal_targets.fa 2> /dev/null && "
   f"kma_index -i acquired_genes.fa -o acquired_genes > /dev/null 2>&1")
log(f"References ready: {len(chrom_targets)} chromosomal genes, {len(quinolone_db) + len(aac6_ib)} acquired-gene alleles")

# ---------------------------------------------------------------------------
# Feature extraction (same functions as notebook 02)
# ---------------------------------------------------------------------------

def read_length_of(url):
    return len(sh(f"set +o pipefail; curl -sL --retry 3 '{url}' | gzip -dc 2>/dev/null "
                  f"| head -n 2 | tail -n 1").strip())


def stream_read_pair(url_r1, url_r2, r1, r2, n_pairs, attempts=4):
    # Stream only the first n_pairs reads of both files; retry failed ENA downloads after a pause,
    # then keep only whole, matching records
    lines = 4 * n_pairs
    fetch = "curl -sL --max-time 1800 '{}' | gzip -dc 2>/dev/null | head -n {} > '{}'"
    for attempt in range(1, attempts + 1):
        sh(f"set +o pipefail; ({fetch.format(url_r1, lines, r1)}) & "
           f"({fetch.format(url_r2, lines, r2)}) & wait")
        complete = min(int(sh(f"wc -l < '{r1}'")), int(sh(f"wc -l < '{r2}'"))) // 4
        if complete >= min(n_pairs, 50_000) or attempt == attempts:
            break
        time.sleep(30 * attempt)
    for path in (r1, r2):
        sh(f"head -n {4 * complete} '{path}' > '{path}.tmp' && mv '{path}.tmp' '{path}'")
    return complete


def good_read(read):
    if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
        return False
    aligned = read.query_alignment_length
    return aligned >= 50 and read.get_tag("NM") <= MAX_MISMATCH * aligned


def consensus_sequence(bam_path, contig, length):
    with pysam.AlignmentFile(bam_path) as bam:
        counts = np.array(bam.count_coverage(
            contig, 0, length, quality_threshold=MIN_BASE_QUAL, read_callback=good_read))
    depth = counts.sum(axis=0)
    majority = counts.max(axis=0)
    agreement = np.divide(majority, depth, out=np.zeros(length), where=depth > 0)
    bases = np.array(list("ACGT"))[counts.argmax(axis=0)]
    bases[(depth < MIN_DEPTH) | (agreement < MIN_BASE_FRAC)] = "N"
    return "".join(bases), depth


def amino_acid_changes(gene, consensus_cds, reference_cds):
    ref_protein = str(Seq(reference_cds).translate(table=11))
    sample_protein = str(Seq(consensus_cds).translate(table=11))
    changes = []
    for position, (ref_aa, sample_aa) in enumerate(zip(ref_protein, sample_protein), start=1):
        if sample_aa not in ("X", ref_aa):
            changes.append(f"{gene}_{ref_aa}{position}{sample_aa}".replace("*", "stop"))
    return changes


def chromosomal_features(r1, r2, work_dir):
    bam = work_dir / "chromosomal.bam"
    sh(f"bwa mem -t {TOOL_THREADS} -v 1 {REF}/chromosomal_targets.fa '{r1}' '{r2}' 2> /dev/null "
       f"| samtools sort -m 256M -o '{bam}' - && samtools index '{bam}'")
    genes = {}
    for gene, info in chrom_targets.items():
        consensus, depth = consensus_sequence(str(bam), gene, info["contig_length"])
        cds_consensus = consensus[FLANK:FLANK + len(info["cds"])]
        cds_depth = depth[FLANK:FLANK + len(info["cds"])]
        genes[gene] = {
            "breadth": round(float((cds_depth >= MIN_DEPTH).mean()), 4),
            "mean_depth": round(float(cds_depth.mean()), 1),
            "changes": amino_acid_changes(gene, cds_consensus, info["cds"]),
        }
    return genes


def gene_family(template):
    name = template.split()[0]
    if "Ib-cr" in name:
        return "aac(6')-Ib-cr"
    if name.startswith("aac(6')"):
        return None  # wild-type aac(6')-Ib: aminoglycoside gene, no quinolone activity
    match = re.match(r"(qnr[A-Z]+|qepA|[Oo]qx[AB])", name)
    return match.group(1).replace("Oqx", "oqx") if match else name.split("_")[0]


def acquired_features(r1, r2, work_dir):
    prefix = work_dir / "kma"
    sh(f"kma -ipe '{r1}' '{r2}' -o '{prefix}' -t_db {REF}/acquired_genes -1t1 -t {TOOL_THREADS} > /dev/null")
    res = pd.read_csv(f"{prefix}.res", sep="\t")
    res.columns = [c.strip().lstrip("#") for c in res.columns]
    hits = []
    for _, row in res.iterrows():
        family = gene_family(row["Template"])
        if family and row["Template_Identity"] >= MIN_IDENTITY and row["Template_Coverage"] >= MIN_COVERAGE:
            hits.append({"family": family, "template": row["Template"].split()[0],
                         "identity": float(row["Template_Identity"]),
                         "coverage": float(row["Template_Coverage"]), "depth": float(row["Depth"])})
    return hits


def mash_sketch(r1, r2, sample_id, work_dir):
    combined = work_dir / "combined.fq"
    sh(f"cat '{r1}' '{r2}' > '{combined}'")
    sh(f"mash sketch -r -m 2 -k 21 -s 10000 -I '{sample_id}' -o '{SKETCH_DIR / sample_id}' '{combined}'")


def process_isolate(row):
    start_time = time.time()
    work_dir = WORK / row.sample_id
    work_dir.mkdir(parents=True, exist_ok=True)
    r1, r2 = work_dir / "r1.fq", work_dir / "r2.fq"
    timings = {}
    step_start = time.time()

    def finish_step(name):
        nonlocal step_start
        timings[name] = round(time.time() - step_start)
        step_start = time.time()

    try:
        read_length = read_length_of(row.url_r1)
        wanted_pairs = int(TARGET_COVERAGE * GENOME_SIZE / (2 * max(read_length, 1)))
        wanted_pairs = min(max(wanted_pairs, 100_000), MAX_READ_PAIRS)
        n_pairs = stream_read_pair(row.url_r1, row.url_r2, r1, r2, wanted_pairs)
        if n_pairs < 50_000:
            raise RuntimeError(f"Too few read pairs downloaded from ENA ({n_pairs}) after 4 attempts")
        if read_length == 0:
            read_length = len(sh(f"sed -n 2p '{r1}'").strip())
        finish_step("download")

        result = {
            "sample_alias": row.sample_alias, "sample_id": row.sample_id,
            "run_accession": row.run_accession, "status": "ok",
            "n_read_pairs": n_pairs, "read_length": read_length,
            "estimated_coverage": round(2 * n_pairs * read_length / GENOME_SIZE, 1),
        }
        result["genes"] = chromosomal_features(r1, r2, work_dir)
        finish_step("mapping")
        result["acquired"] = acquired_features(r1, r2, work_dir)
        finish_step("kma")
        mash_sketch(r1, r2, row.sample_id, work_dir)
        finish_step("mash")
    except Exception as error:
        result = {"sample_alias": row.sample_alias, "sample_id": row.sample_id,
                  "status": "error", "message": str(error)[-1500:]}
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    result["timings"] = timings
    result["seconds"] = round(time.time() - start_time, 1)
    tmp_path = RESULTS_DIR / f"{row.sample_id}.json.tmp"
    tmp_path.write_text(json.dumps(result, indent=1))
    tmp_path.replace(RESULTS_DIR / f"{row.sample_id}.json")
    return result


def describe(result):
    steps = ", ".join(f"{k} {v}s" for k, v in result["timings"].items())
    if result["status"] == "ok":
        mutations = [c for g in result["genes"].values() for c in g["changes"]
                     if c.split("_")[0] in ("gyrA", "parC", "parE")]
        acquired = sorted({h["family"] for h in result["acquired"]})
        return (f"[{steps}] cov {result['estimated_coverage']}x | "
                f"gyrA/parC/parE: {', '.join(mutations) or '-'} | acquired: {', '.join(acquired) or '-'}")
    lines = [l for l in result["message"].splitlines() if l.strip()]
    return f"[{steps}] ERROR: " + " / ".join(lines[-3:])[-400:]


def already_done(sample_id):
    path = RESULTS_DIR / f"{sample_id}.json"
    if not path.exists():
        return False
    return json.loads(path.read_text())["status"] == "ok" or args.no_retry_errors


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if not args.combine_only:
    to_process = isolates if args.n_isolates is None else isolates.head(args.n_isolates)
    to_process = to_process[to_process["url_r1"].notna()]
    pending = [row for row in to_process.itertuples() if not already_done(row.sample_id)]
    log(f"{len(to_process) - len(pending)} already done, {len(pending)} to process, {N_WORKERS} in parallel")

    run_start = time.time()
    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(process_isolate, row): row for row in pending}
        for i, future in enumerate(as_completed(futures), start=1):
            row, result = futures[future], future.result()
            elapsed = (time.time() - run_start) / 60
            remaining = elapsed / i * (len(pending) - i)
            log(f"[{i}/{len(pending)}] {row.sample_alias} ({result['seconds']:.0f}s, "
                f"{elapsed:.0f} min elapsed, ~{remaining:.0f} min left) {describe(result)}")

# ---------------------------------------------------------------------------
# Combine into feature tables
# ---------------------------------------------------------------------------
results = [json.loads(p.read_text()) for p in sorted(RESULTS_DIR.glob("*.json"))]
ok_results = [r for r in results if r["status"] == "ok" and (SKETCH_DIR / f"{r['sample_id']}.msh").exists()]
failed = [r for r in results if r["status"] == "error"]
log(f"Successful: {len(ok_results)}, failed: {len(failed)}")
for r in failed:
    log(f"   {r['sample_alias']} -> {r['message'].splitlines()[-1][:150]}")

feature_rows, qc_rows = [], []
for r in ok_results:
    features = {"sample_alias": r["sample_alias"]}
    for gene in r["genes"].values():
        for change in gene["changes"]:
            features[change] = 1
    for hit in r["acquired"]:
        features[f"acq_{hit['family']}"] = 1
    feature_rows.append(features)

    qc = {k: r[k] for k in ["sample_alias", "run_accession", "n_read_pairs",
                            "read_length", "estimated_coverage", "seconds"]}
    for gene_name, gene in r["genes"].items():
        qc[f"{gene_name}_breadth"] = gene["breadth"]
        qc[f"{gene_name}_depth"] = gene["mean_depth"]
    qc["acquired_alleles"] = ";".join(h["template"] for h in r["acquired"])
    qc_rows.append(qc)

features = pd.DataFrame(feature_rows).fillna(0)
feature_cols = sorted(c for c in features.columns if c != "sample_alias")
features[feature_cols] = features[feature_cols].astype(int)
features[["sample_alias"] + feature_cols].to_csv(OUT / "features.csv", index=False)
pd.DataFrame(qc_rows).to_csv(OUT / "qc.csv", index=False)

ok_ids = [r["sample_id"] for r in ok_results]
(WORK / "sketch_list.txt").write_text("\n".join(str(SKETCH_DIR / f"{i}.msh") for i in ok_ids))
sh(f"rm -f {WORK}/all_sketches.msh && mash paste {WORK}/all_sketches -l {WORK}/sketch_list.txt")
dist_text = sh(f"mash dist -p {os.cpu_count()} -t {WORK}/all_sketches.msh {WORK}/all_sketches.msh")
distances = pd.read_csv(io.StringIO(dist_text), sep="\t", index_col=0)
distances.columns = distances.index
id_to_alias = dict(zip(isolates["sample_id"], isolates["sample_alias"]))
distances = distances.rename(index=id_to_alias, columns=id_to_alias)
distances.index.name = "sample_alias"
distances.to_csv(OUT / "mash_distances.csv")

log(f"Written to {OUT}: features.csv {features.shape}, qc.csv, mash_distances.csv {distances.shape}")
