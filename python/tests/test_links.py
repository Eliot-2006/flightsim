"""Every relative link and image in the published Markdown documents must resolve to a file in the repository."""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DOCS = ["README.md", "PROGRESS.md", "docs/PRD.md", "docs/design.md", "docs/results.md", "data/f16/SOURCES.md"]
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")


@pytest.mark.parametrize("doc", DOCS)
def test_relative_links_resolve(doc):
    path = ROOT / doc
    broken = []
    for target in LINK.findall(path.read_text()):
        if re.match(r"[a-z]+:", target) or target.startswith("#"):
            continue  # external URLs are checked by hand at release time; anchors are page-internal
        file = target.split("#", 1)[0]
        if not (path.parent / file).resolve().exists():
            broken.append(target)
    assert not broken, f"{doc}: broken relative links {broken}"
