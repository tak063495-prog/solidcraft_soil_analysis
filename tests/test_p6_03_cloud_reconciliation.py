"""Include P6-03 reconciliation in the ordinary Python inventory."""
import importlib.util
from pathlib import Path
_PATH=Path(__file__).resolve().parents[1]/'interim/2026-12-07-p6-03-cloud-reconciliation/test_audit.py'
_SPEC=importlib.util.spec_from_file_location('p603_reconciliation_cases',_PATH)
_CASES=importlib.util.module_from_spec(_SPEC);_SPEC.loader.exec_module(_CASES)
P603AuditTests=_CASES.P603AuditTests
