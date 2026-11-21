"""Verify SHA256, extract to a new folder, and remeasure without any solver."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from study import measure,write,serializable
HERE=Path(__file__).resolve().parent

def verify(output):
    manifest=json.loads((HERE/'archive-sha256.json').read_text());archive=HERE/'evidence.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==manifest['archive_sha256']
    output.mkdir(parents=True,exist_ok=False)
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(manifest['entries'])
        for n,h in manifest['entries'].items():
            raw=z.read(n);assert hashlib.sha256(raw).hexdigest()==h,('changed entry',n)
            p=output/n;assert p.resolve().is_relative_to(output.resolve());p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
    actual=serializable(measure(output/'raw/paths'));expected=json.loads((HERE/'measurements.json').read_text())
    assert actual==expected, 'measurements differ'
    report=dict(software_check_passed=True,all_entries_sha256_verified=len(manifest['entries']),all_measurements_exactly_reproduced=True,
                solver_executed=False,passed=False,history_converged=False,windows_ci='unverified')
    write(output/'remeasurement.json',report);print(json.dumps(report));return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.output.resolve())
