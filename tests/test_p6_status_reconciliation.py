"""Include the independent P6 ledger audit in the ordinary Python inventory."""
import importlib.util
from pathlib import Path

_PATH = (Path(__file__).resolve().parents[1]
         / 'interim/2026-12-06-p6-01-status-reconciliation/test_status_audit.py')
_SPEC = importlib.util.spec_from_file_location('p6_status_reconciliation_cases', _PATH)
_CASES = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_CASES)
AuditMutationTests = _CASES.AuditMutationTests
FixedRepositoryTests = _CASES.FixedRepositoryTests
