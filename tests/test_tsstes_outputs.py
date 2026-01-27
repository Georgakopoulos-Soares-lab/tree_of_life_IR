import gzip
from pathlib import Path
import pandas as pd
import pytest
import numpy as np

DENSITY_DIRS = [
    Path("nonbdna_pipeline/extractions_IR/tss_tes_density"),
    Path("extractions_IR/tss_tes_density"),
]
DENSITY_FILE_PATTERN = "tss_tes_density_IR_bucket_0.tsv.gz"
META_COLS = {"#assembly_accession", "pattern", "site", "partition", "biotype", "polarity"}

def find_density_file() -> Path:
    for d in DENSITY_DIRS:
        p = d / DENSITY_FILE_PATTERN
        if p.is_file():
            return p
    pytest.skip(f"Density file {DENSITY_FILE_PATTERN} not found in known dirs: {DENSITY_DIRS}")

def load_density_df(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_table(f)
    return df

@pytest.fixture 
def window_size() -> int:
    return 10

@pytest.fixture
def expected_vectors_dict(window_size: int) -> dict[str, dict[int, int]]:
    """
    Quiz answers: expected per-position vectors for TSS/TES with W=10 based on user's notes.
    Keys: "TSS" and "TES". Positions are integers in [-W, +W].
    """
    # Expected vectors built from strand-aware per-chromosome overlaps:
    # NC_008229.1:
    #   - gene 1: 1..35 (-), TSS at 35 → overlaps IR (10..25) only at p=25 → offset +10 → index -1 gets +1
    #   - gene 2: 25..55 (+), TSS at 25:
    #       IR (10..25) contributes p=15..25 → offsets -10..0
    #       IR (5..20)  contributes p=15..20 → offsets -10..-5
    #     net TSS counts modeled as: indices 0..4 (+1 each), indices 5..10 (+2 each)
    # NC_008229.2: no overlaps
    # NC_008229.3: 10..11 (-), TSS at 11 → IR (10..16) contributes 6 bp above TSS → last 6 bins get +1
    answers = dict()
    expected_vector_TSS = np.zeros(window_size * 2 + 1, dtype=int)
    expected_vector_TSS[:5] += 1
    expected_vector_TSS[5:11] += 2
    expected_vector_TSS[-1] += 1

    expected_vector_TES = np.zeros(window_size * 2 + 1, dtype=int) 
    expected_vector_TES[:10] += 1
    expected_vector_TES[4:6] += 1
    answers.update({"GCF_test.1": {"TSS": list(expected_vector_TSS), "TES": list(expected_vector_TES)}})

    expected_vector_TSS = np.zeros(window_size * 2 + 1, dtype=int)
    expected_vector_TES = np.zeros(window_size * 2 + 1, dtype=int)
    answers.update({"GCF_test.2": {"TSS": list(expected_vector_TSS), "TES": list(expected_vector_TES)}})
    return answers


@pytest.mark.parametrize("accession_prefix, expected_tss_sum, expected_tes_sum", [
    ("GCF_test.1", 18, 12),
    ("GCF_test.2", 0, 0),
])
def test_accession_vector_sums(accession_prefix, expected_tss_sum, expected_tes_sum):
    path = find_density_file()
    df = load_density_df(path)
    window_size = 10
    cols = list(map(str, range(-window_size, window_size+1)))

    rows = df[df["#assembly_accession"].astype(str).str.contains(accession_prefix)]
    if rows.empty:
        pytest.skip(f"No rows found for accession prefix {accession_prefix} in {path}")

    tss = rows[rows["site"] == "TSS"]
    tes = rows[rows["site"] == "TES"]
    assert not tss.empty, f"Missing TSS row for {accession_prefix}"
    assert not tes.empty, f"Missing TES row for {accession_prefix}"

    tss_sum = int(tss[cols].to_numpy().sum())
    tes_sum = int(tes[cols].to_numpy().sum())

    assert tss_sum == expected_tss_sum, f"TSS sum mismatch for {accession_prefix}: got {tss_sum}, expected {expected_tss_sum}"
    assert tes_sum == expected_tes_sum, f"TES sum mismatch for {accession_prefix}: got {tes_sum}, expected {expected_tes_sum}"

@pytest.mark.parametrize("accession_prefix", ["GCF_test.1", "GCF_test.2"])
def test_accession_vector_identity(accession_prefix, expected_vectors_dict):
    path = find_density_file()
    df = load_density_df(path)
    window_size = 10
    cols = list(map(str, range(-window_size, window_size+1)))

    rows = df[df["#assembly_accession"].astype(str).str.contains(accession_prefix)]
    if rows.empty:
        pytest.skip(f"No rows found for accession prefix {accession_prefix} in {path}")

    tss = rows[rows["site"] == "TSS"]
    tes = rows[rows["site"] == "TES"]
    assert not tss.empty, f"Missing TSS row for {accession_prefix}"
    assert not tes.empty, f"Missing TES row for {accession_prefix}"

    expected_tss_array = expected_vectors_dict[accession_prefix]["TSS"]
    expected_tes_array = expected_vectors_dict[accession_prefix]["TES"]

    assert expected_tss_array == list(tss[cols].to_numpy().flatten()), f"TSS vector mismatch for {accession_prefix}"
    assert expected_tes_array == list(tes[cols].to_numpy().flatten()), f"TES vector mismatch for {accession_prefix}"


def test_print_quiz_answers(expected_vectors_dict):
    """Emit the expected vectors as dictionaries for quick inspection."""
    assert sum(expected_vectors_dict["GCF_test.1"]["TSS"]) == 18
    assert sum(expected_vectors_dict["GCF_test.1"]["TES"]) == 12
    assert sum(expected_vectors_dict["GCF_test.2"]["TSS"]) == 0
    assert sum(expected_vectors_dict["GCF_test.2"]["TES"]) == 0