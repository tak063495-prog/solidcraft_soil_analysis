"""Read the exact production blob used by a historical diagnosis."""
import hashlib
import subprocess


def diagnostic_bytes(root, commit, path):
    # Keep checking retained reports, workflow rules and numerical evidence in the checkout.
    if not path.startswith(('src/', 'examples/', 'tools/')):
        return (root / path).read_bytes()

    def git(*args):
        return subprocess.check_output(['git', '-C', str(root), *args])

    blob = git('rev-parse', '--verify', f'{commit}:{path}').decode('ascii').strip()
    data = git('cat-file', 'blob', blob)
    actual = hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()
    if actual != blob:
        raise ValueError(f'{commit}:{path} blob mismatch: {actual} != {blob}')
    return data
