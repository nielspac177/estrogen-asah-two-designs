"""P0.4 — config loader loads the shipped example template."""

from estrogen_asah_dci import __version__
from estrogen_asah_dci.config import REPO_ROOT, DataPaths, load_paths


def test_package_imports():
    assert __version__ == "0.1.0"


def test_load_example_paths():
    paths = load_paths(REPO_ROOT / "config" / "paths.example.yaml")
    assert isinstance(paths, DataPaths)
    # template placeholder values, not real machine paths
    assert str(paths.mimic_iv_db).endswith("mimic4.db")
    assert "eicu" in str(paths.eicu_dir).lower()


def test_example_paths_do_not_exist():
    # the template must not accidentally resolve to real data
    paths = load_paths(REPO_ROOT / "config" / "paths.example.yaml")
    assert paths.exists() == {"mimic_iv_db": False, "eicu_dir": False}
