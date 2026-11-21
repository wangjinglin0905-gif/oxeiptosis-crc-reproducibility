"""Audit and prepare GSE294385 counts and H&E labels without computing scores."""

import collections
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
import time

import numpy as np
import openpyxl
import pandas as pd
import requests
from scipy.io import mmread
from threadpoolctl import threadpool_limits

import os
ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('OX4_OUTPUT', ROOT/'run'))
GEO = Path(os.environ.get('OX4_INPUTS', ROOT/'inputs'))
LABEL_DIR = GEO
PIN = 'b6c40903eeda67d812712f2bb77064041db38e6c'
ANNOTATIONS = "visium_meta_after_qc.rds"
LABEL_URL = f"https://raw.githubusercontent.com/yliuup/CRC_micromets_ST/{PIN}/Meta_data/{ANNOTATIONS}"
S1_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S1535610826002965-mmc2.xlsx"
GEO_FILES = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294385/suppl/filelist.txt"
SAMPLES = {
    "M-ST-13": (1, "GSM8903668"), "M-ST-15": (6, "GSM8903670"),
    "M-ST-29": (11, "GSM8903684"), "M-ST-30": (11, "GSM8903685"),
    "M-ST-31": (9, "GSM8903686"), "M-ST-32": (9, "GSM8903687"),
    "M-ST-33": (10, "GSM8903688"), "M-ST-34": (10, "GSM8903689"),
}
SETS = {
    "OX4": ("KEAP1", "PGAM5", "AIFM1", "OTUD1"),
    "epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19"),
    "stromal": ("COL1A1", "COL1A2", "COL3A1", "DCN"),
}
ZONES = ("mucosa and submucosa", "muscularis propria", "subserosa", "invasive front/border")
TUMOR = dict(zip([f"Tumor cells in {z}" for z in ZONES[:3]]
                 + ["Tumor cells at tumor invasive front/border"], ZONES))
STROMA = dict(zip([f"Tumor stroma in {z}" for z in ZONES[:3]]
                  + ["Tumor stroma at tumor invasive front/border"], ZONES))
SEED = 20261118


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def inputs():
    """Verify sample selection against publisher table, GEO catalogue and download receipts."""
    counts = pd.read_csv(io.StringIO(request(GEO_FILES).decode()), sep="\t")
    file_sizes = counts.set_index("Name")["Size"].to_dict()
    book_bytes = request(S1_URL)
    book = openpyxl.load_workbook(io.BytesIO(book_bytes), read_only=True, data_only=True)
    patient = tissue = None
    primary = {}
    for row_number, row in enumerate(book["Tab 1"].values, 1):
        if row_number < 4:
            continue
        if row[0] is not None:
            patient, tissue = row[0], None
        if row[2] is not None:
            tissue = str(row[2]).strip()
        if tissue == "Primary CRC" and str(row[1]).startswith("M-ST-") and row[3] == "✓":
            primary[row[1]] = {"patient": int(patient), "row": row_number}
    assert {key: value["patient"] for key, value in primary.items()} == {
        key: value[0] for key, value in SAMPLES.items()
    }
    label_hash = sha(LABEL_DIR / ANNOTATIONS)
    assert label_hash == hashlib.sha256(request(LABEL_URL)).hexdigest()
    manifest = json.loads((ROOT/"tables/S40_input_provenance.json").read_text(encoding="utf8"))
    assert sha(LABEL_DIR / ANNOTATIONS)==manifest["annotation_sha256"]
    assert hashlib.sha256(book_bytes).hexdigest()==manifest["supplement_sha256"]
    rows = []
    for sample, (patient, gsm) in SAMPLES.items():
        name = f"{gsm}_{sample}.tar.gz"
        archive = GEO / name
        expected = next(r for r in manifest['archives'] if r['sample'] == sample)
        digest = sha(archive)
        assert archive.stat().st_size == file_sizes[name] == expected['size_bytes']
        assert digest == expected['sha256'] and gsm == expected['gsm']
        rows.append(dict(sample=sample, patient=patient, gsm=gsm, path=str(archive),
                         size_bytes=archive.stat().st_size, sha256=digest, url=expected['url'],
                         table_s1_row=primary[sample]['row']))
    return dict(supplement_url=S1_URL, supplement_sha256=hashlib.sha256(book_bytes).hexdigest(),
                source_filelist=GEO_FILES, annotation_url=LABEL_URL,
                annotation_sha256=label_hash, archives=rows)


def request(url):
    response = requests.get(url, timeout=35)
    response.raise_for_status()
    return response.content


def matrix_archive(path, sample):
    root = f"{sample}/"
    names = ["filtered_feature_bc_matrix/" + name for name in
             ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")]
    names += ["spatial/tissue_positions.csv", "spatial/scalefactors_json.json",
              "spatial/tissue_hires_image.png"]
    files = {}
    with tarfile.open(path, "r|gz") as bundle:
        for member in bundle:
            suffix = member.name.removeprefix(root)
            if member.isfile() and member.name.startswith(root) and suffix in names:
                assert suffix not in files, member.name
                files[suffix] = bundle.extractfile(member).read()
    assert set(files) == set(names), (sample, set(names) - set(files))
    barcodes = gzip.decompress(files[names[0]]).decode().splitlines()
    features = pd.read_csv(io.BytesIO(gzip.decompress(files[names[1]])), sep="\t", header=None)
    with gzip.GzipFile(fileobj=io.BytesIO(files[names[2]])) as matrix:
        with threadpool_limits(limits=4):
            blocks = mmread(matrix)
        assert not matrix.read(1), f"Unused matrix content in {sample}"
    assert blocks.shape == (len(features), len(barcodes)), sample
    assert np.issubdtype(blocks.data.dtype, np.integer) and (blocks.data > 0).all(), sample
    X = blocks.tocsr()
    assert X.nnz == blocks.nnz, f"Duplicate matrix coordinates in {sample}"
    del blocks
    assert features[0].is_unique and len(set(barcodes)) == len(barcodes)
    positions = pd.read_csv(io.BytesIO(files[names[3]])).set_index("barcode", verify_integrity=True)
    positions = positions.loc[barcodes].copy()
    scale = json.loads(files[names[4]])["tissue_hires_scalef"]
    # PNG IHDR dimensions, kept as a registration bounds check without re-encoding the image.
    image = files[names[5]]
    assert image[:8] == b"\x89PNG\r\n\x1a\n" and image[12:16] == b"IHDR"
    width, height = [int.from_bytes(image[n:n+4], "big") for n in (16, 20)]
    positions["x_hires"] = positions.pxl_col_in_fullres * scale
    positions["y_hires"] = positions.pxl_row_in_fullres * scale
    return X, features, positions, width, height


def background_constructs(symbols, X, library_sizes, sd, available, sample):
    """Use the frozen expression/detection matching rule in supplementary Text S4."""
    assert X.shape[1] > 1
    mean = np.asarray(X.multiply(10000 / library_sizes).mean(axis=1)).ravel()
    detect = np.diff(X.indptr) / X.shape[1]
    occurrences = collections.Counter(symbols)
    exclude = set().union(*SETS.values())
    pool = np.array([i for i, gene in enumerate(symbols)
                     if occurrences[gene] == 1 and gene not in exclude
                     and not gene.startswith(("MT-", "RPL", "RPS")) and sd[i] > 1e-12])
    sets = available.copy()
    matches, balance, pool_sizes = [], [], {}
    for gene in SETS["OX4"]:
        idx = symbols.index(gene)
        keep = ((mean[pool] >= .5 * mean[idx]) & (mean[pool] <= 2 * mean[idx])
                & (abs(detect[pool] - detect[idx]) <= .05))
        valid = pool[keep]
        pool_sizes[gene] = len(valid)
        distance = abs(np.log2(mean[valid] / mean[idx])) + abs(detect[valid] - detect[idx]) / .05
        nearest = valid[np.lexsort((valid, distance))[:30]]
        matches.append(nearest)
        balance.extend(dict(target=gene, candidate=symbols[i], mean_ratio=mean[i] / mean[idx],
                            detection_difference=detect[i] - detect[idx]) for i in nearest)
    if min(pool_sizes.values()) >= 10:
        rng = np.random.default_rng(np.random.SeedSequence([SEED, list(SAMPLES).index(sample), 99]))
        seen = set()
        for _ in range(100000):
            chosen = tuple(int(rng.choice(pool)) for pool in matches)
            key = tuple(sorted(chosen))
            if len(set(chosen)) == 4 and key not in seen:
                seen.add(key)
                sets[f"random_{len(seen):03d}"] = list(chosen)
                if len(seen) == 100:
                    break
        if len(seen) != 100:
            return available.copy(), balance, pool_sizes
    return sets, balance, pool_sizes


def prepare(source, labels):
    start = time.monotonic()
    sample = source["sample"]
    destination = OUT / "sections" / sample
    destination.mkdir(parents=True, exist_ok=True)
    X, features, positions, width, height = matrix_archive(Path(source["path"]), sample)
    available_labels = labels[labels["sample"].eq(sample)].set_index("spot", verify_integrity=True)
    assert len(available_labels) and available_labels.index.isin(positions.index).all()
    data = positions.join(available_labels.drop(columns="sample"), validate="one_to_one")
    data["has_label"] = data["Layer1"].eq("Primary Colon")
    data["inside_image"] = data.x_hires.between(0, width, inclusive="left") & data.y_hires.between(
        0, height, inclusive="left")
    types = features[2].eq("Gene Expression").to_numpy()
    symbols = features.loc[types, 1].tolist()
    X = X[types].tocsr()
    total = np.asarray(X.sum(axis=0)).ravel()
    data["library"] = total
    any_target = data.Layer3.isin(set(TUMOR) | set(STROMA))
    primary = data["Layer2"].eq("Primary colon carcinoma")
    data["class"] = np.select([primary & data.Layer3.isin(TUMOR), primary & data.Layer3.isin(STROMA)],
                              ["tumor", "stroma"], default="other")
    data["zone"] = data.Layer3.map(TUMOR | STROMA)
    usable = data.has_label & data.in_tissue.eq(1) & data.library.ge(500) & data.inside_image
    usable &= data.Layer2.notna() & data.Layer3.notna()
    usable &= ~data.Layer2.isin(["Other", ""]) & ~data.Layer3.isin(["Other", ""])
    data["eligible"] = usable
    data.to_csv(destination / "spot_eligibility.csv.gz", compression="gzip")
    layer_counts = data.groupby(["Layer1", "Layer2", "Layer3", "class", "eligible"],
                                dropna=False).size().rename("n").reset_index()
    layer_counts.to_csv(destination / "label_counts.csv", index=False)
    number = data.loc[usable, "class"].value_counts()

    X = X[:, usable.to_numpy()].tocsr()
    library_sizes = total[usable]
    normalized = X.multiply(10000 / library_sizes).tocsr()
    normalized.data = np.log1p(normalized.data)
    center = np.asarray(normalized.mean(axis=1)).ravel()
    variance = np.asarray(normalized.power(2).mean(axis=1)).ravel() - center**2
    sd = np.sqrt(np.maximum(variance, 0) * (X.shape[1] / (X.shape[1] - 1)))
    del normalized
    coverage, constructs = [], {}
    for label, genes in SETS.items():
        indices = []
        for gene in genes:
            hits = [i for i, symbol in enumerate(symbols) if symbol == gene]
            index = hits[0] if len(hits) == 1 else None
            passes = index is not None and sd[index] > 1e-12
            coverage.append(dict(construct=label, gene=gene, matrix_features=len(hits),
                                 feature_id="" if index is None else features.loc[types, 0].iloc[index],
                                 detected_spots=0 if index is None else int(X[index].nnz),
                                 variation=bool(passes)))
            if passes:
                indices.append(index)
        if len(indices) == len(genes):
            constructs[label] = indices
    pd.DataFrame(coverage).to_csv(destination / "gene_coverage.csv", index=False)
    gate = bool("OX4" in constructs and min(number.get("tumor", 0), number.get("stroma", 0)) >= 30)
    pool_sizes = {}
    if gate:
        constructs, balance, pool_sizes = background_constructs(symbols, X, library_sizes, sd,
                                                                 constructs, sample)
        pd.DataFrame(balance).to_csv(destination / "match_candidates.csv", index=False)
        union = sorted(set().union(*constructs.values()))
        column = {index: i for i, index in enumerate(union)}
        weights = np.zeros((len(union), len(constructs)))
        for i, indices in enumerate(constructs.values()):
            weights[[column[index] for index in indices], i] = .25
        pd.DataFrame([dict(construct=label, **{f"gene{i}": symbols[index]
                                              for i, index in enumerate(indices, 1)})
                      for label, indices in constructs.items()]).to_csv(
                          destination / "constructs.csv", index=False)
        np.savez_compressed(destination / "eligible_counts.npz", count=X[union].toarray().T,
                            library=library_sizes, center=center[union], sd=sd[union],
                            weights=weights, names=np.asarray(list(constructs)),
                            compartments=data.loc[usable, "class"].to_numpy(dtype=str),
                            zones=data.loc[usable, "zone"].fillna("").to_numpy(dtype=str),
                            barcodes=data.index[usable].to_numpy(dtype=str))

    result = dict(sample=sample, patient=source["patient"], raw_spots=len(data),
                  source_labels=len(available_labels), raw_spots_outside_source_labels=int(
                      (~data.has_label).sum()), eligible_spots=int(usable.sum()),
                  labelled_outside_hires_image=int((data.has_label & ~data.inside_image).sum()),
                  labelled_outside_hires_image_at_ge500=int(
                      (data.has_label & ~data.inside_image & data.library.ge(500)).sum()),
                  n_tumor=int(number.get("tumor", 0)), n_stroma=int(number.get("stroma", 0)),
                  target_layer_discrepancies=int((any_target & ~primary).sum()),
                  target_layer_discrepancies_among_eligible=int((any_target & ~primary & usable).sum()),
                  gene_expression_features=len(symbols), nonzero_entries=int(X.nnz),
                  usable_in_tissue_ge500_before_other_filter=int(
                      (data.has_label & data.in_tissue.eq(1) & data.library.ge(500)).sum()),
                  position_resolution=[width, height], complete_constructs=list(constructs)[:3],
                  background_count=sum(key.startswith("random_") for key in constructs),
                  match_pool_sizes=pool_sizes, comparable=gate, prepare_seconds=time.monotonic()-start)
    write_json(destination / "gate.json", result)
    print(sample, result["eligible_spots"], result["n_tumor"], result["n_stroma"],
          result["background_count"], round(result["prepare_seconds"], 1), flush=True)
    return result


def main():
    provenance = inputs()
    write_json(OUT / "input_provenance.json", provenance)
    labels = pd.read_csv(OUT / "source_primary_labels.csv", dtype=str)
    assert len(labels) == 60801 and labels.Barcode.is_unique
    labels[["sample", "spot"]] = labels.Barcode.str.split("~", n=1, expand=True)
    assert sorted(labels["sample"].unique()) == list(SAMPLES)
    gates = []
    for record in provenance["archives"]:
        try:
            gates.append(prepare(record, labels))
        except Exception as error:
            gates.append(dict(sample=record["sample"], patient=record["patient"],
                              comparable=False, error=f"{type(error).__name__}: {error}"))
            print(gates[-1], flush=True)
        write_json(OUT / "input_gate.json", gates)
    assert sum(row["comparable"] for row in gates) > 0, "No eligible comparison"


if __name__ == "__main__":
    main()
