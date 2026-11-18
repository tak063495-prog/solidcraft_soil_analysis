"""Check the working bytes against pre-18d54b1 preservation manifests."""
import hashlib


WORKFLOW = '.github/workflows/fem-ci.yml'
# Pin the complete approved line, including its indentation and LF.
PC_EVIDENCE_RULE = (
    b"            ${{ (matrix.group == 'python-core-0' || matrix.group == 'python-core-1')"
    b" && '/interim/2027-01-08-pc-failure-diagnosis/evidence.zip.part*'"
    b" || '!/interim/2027-01-08-pc-failure-diagnosis/evidence.zip.part*' }}\n"
)
NEXT_STEP = b"      - name: Set up Python\n"


def prior_record_sha256(name, data):
    if name == WORKFLOW:
        lines = data.splitlines(keepends=True)
        if lines.count(PC_EVIDENCE_RULE) != 1:
            raise ValueError('expected exactly one approved PC evidence sparse-checkout line')
        index = lines.index(PC_EVIDENCE_RULE)
        if lines[index + 1:index + 2] != [NEXT_STEP]:
            raise ValueError('PC evidence rule moved from the end of sparse-checkout')
        del lines[index]
        data = b''.join(lines)
    return hashlib.sha256(data).hexdigest()
