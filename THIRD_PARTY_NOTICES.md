# Third-party references

SolidCraft integration changes target https://github.com/tangency-labs/solidcraft at the commit recorded in `sources.lock.json`. Its MIT/Apache-2.0 notices, attribution and third-party license inventory accompany the patch in `integrations/solidcraft/`.

The 2D comparison target is https://github.com/tak063495-prog/excel-vba-soil-fem at the pinned revision in `sources.lock.json`. Its full source and Excel installation are not bundled. The pinned Griffiths 661-node coordinates, supports and published displacement values are redistributed in tests/data/griffiths_2d_reference.json under Apache-2.0; see tests/data/GRIFFITHS_REFERENCE_LICENSE.txt and docs/GRIFFITHS_2D_COMPARISON.ja.md for attribution and provenance. Rust and Python dependencies are declared in `Cargo.lock` and `requirements.lock.txt`; dependency runtimes are installed separately.
