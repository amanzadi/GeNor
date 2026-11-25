import json
from pathlib import Path

import pytest

GOLDEN = Path(__file__).parent / "golden"


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """Every test gets its own GeNor data/config/jobs dirs."""
    monkeypatch.setenv("GENOR_HOME", str(tmp_path / "genor_home"))
    monkeypatch.delenv("GENOR_DATA_DIR", raising=False)
    monkeypatch.delenv("GENOR_CONFIG_DIR", raising=False)
    monkeypatch.delenv("GENOR_JOBS_DIR", raising=False)
    yield


@pytest.fixture(scope="session")
def expected_r():
    return json.loads((GOLDEN / "expected_r.json").read_text())


@pytest.fixture
def example_dir(tmp_path):
    from genor.commands import example

    example(tmp_path / "ex")
    return tmp_path / "ex"


@pytest.fixture
def loaded_pair(example_dir):
    from genor.commands import sumstats_load

    sumstats_load(example_dir / "bmi_sim.tsv.gz", "BMI_sim", example_dir / "bmi", ancestry="EUR", plots=False)
    sumstats_load(example_dir / "mdd_sim.tsv.gz", "MDD_sim", example_dir / "mdd", ancestry="EUR",
                  population_prevalence=0.15, plots=False)
    return example_dir / "bmi", example_dir / "mdd"
