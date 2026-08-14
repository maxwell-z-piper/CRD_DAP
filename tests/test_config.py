from pathlib import Path

from crd_utils.config import load_config


def test_template_config_loads_without_strict_paths():
    path = Path(__file__).resolve().parents[1] / "config" / "target_config_template.py"
    cfg = load_config(path, validate=False)
    assert hasattr(cfg, "TARGET_NAME")
    assert hasattr(cfg, "BL_MASTER_ARC")
