"""Exercise the workflow's sparse patterns with real Git, including Windows text I/O."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import run_ci_python_evidence_shards as runner

HERE = ROOT/'interim/2026-12-14-ci-timeout-fa22e14'
FOLDERS = ('2026-12-11-c01-global-step-control',
           '2026-12-12-c01-switch-event-localization',
           '2026-12-13-c01-cause-diagnosis')


class CITimeoutConfiguration(unittest.TestCase):
    def workflow(self):
        return (ROOT/'.github/workflows/fem-ci.yml').read_text(encoding='utf-8')

    def test_workflow_runs_every_python_partition_and_keeps_exact_history(self):
        text = self.workflow()
        groups = re.search(r'        group: \[([^\]]+)\]', text)[1].split(', ')
        self.assertEqual({g for g in groups if g.startswith('python-')}, set(runner.GROUPS))
        self.assertIn('fetch-depth: 0', text)
        self.assertLess(text.index('git config --global core.autocrlf false'),
                        text.index('name: Check out sources'))
        self.assertIn('fail-fast: false', text)
        self.assertIn("$env:GROUP_RESULT -ne 'success'", text)

    def test_real_git_sparse_profiles_keep_all_sources_and_each_required_archive(self):
        text = self.workflow()
        block = text.split('          sparse-checkout: |\n', 1)[1].split('      - name:', 1)[0]
        templates = [line.strip() for line in block.splitlines() if line.strip()]
        self.assertEqual(len(templates), 10)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def git(*args, **kwargs):
                return subprocess.run(['git', '-C', str(root), *args], check=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
            git('init', '--quiet')
            git('config', 'user.name', 'CI fixture')
            git('config', 'user.email', 'ci-fixture@example.invalid')
            sources = ['README.md', 'src/lib.rs', 'tests/test_example.py',
                       'docs/status.ja.md', 'interim/unrelated/validation.json']
            archives = ['interim/'+folder+'/evidence.zip.part000' for folder in FOLDERS]
            unrelated = ['interim/unrelated/evidence.tar.xz.part000',
                         'interim/unrelated/evidence.zip', 'interim/unrelated/history.json.gz',
                         'interim/unrelated/data.npz']
            prototypes = 'interim/2026-12-12-c01-switch-event-localization/retained-failures/prototypes.zip'
            omitted = 'interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part000'
            for name in sources+archives+unrelated+[prototypes, omitted]:
                path = root/name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b'fixture\n')
            git('add', '.')
            git('commit', '--quiet', '-m', 'fixture')
            for group in (*runner.GROUPS, 'rust-core', 'rust-integration-0',
                          'rust-integration-1', 'verification'):
                patterns = []
                for line in templates:
                    if line.startswith('${{'):
                        match = re.fullmatch(r"\$\{\{ \((.+)\) && '([^']+)' \|\| '([^']+)' \}\}", line)
                        self.assertIsNotNone(match)
                        condition, yes, no = match.groups()
                        owners = re.findall(r"matrix.group == '([^']+)'", condition)
                        self.assertEqual(condition, ' || '.join("matrix.group == '"+g+"'" for g in owners))
                        patterns.append(yes if group in owners else no)
                    else:
                        patterns.append(line)
                git('sparse-checkout', 'set', '--no-cone', '--stdin',
                    input=('\n'.join(patterns)+'\n').encode('utf-8'))
                for name in sources:
                    self.assertTrue((root/name).is_file(), (group, name))
                expected = {'python-c01-switch': {0}, 'python-c01-path': {1},
                            'python-c01-cause': {0, 1, 2}}.get(group, set())
                self.assertEqual({i for i, name in enumerate(archives) if (root/name).is_file()}, expected, group)
                for name in unrelated:
                    self.assertEqual((root/name).is_file(), not group.startswith('python-c01-'), (group, name))
                self.assertEqual((root/prototypes).is_file(), group in ('python-c01-path', 'python-c01-cause'))
                self.assertFalse((root/omitted).exists())

    def test_new_document_decisions_survive_real_crlf_bytes(self):
        for name in ('README.ja.md', 'REPRODUCE.ja.md'):
            text = (HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('C-01は未達、passed=false、原因は判定不能。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)/name
                path.write_bytes(text.replace('\n', '\r\n').encode('utf-8'))
                self.assertEqual(path.read_text(encoding='utf-8'), text)


if __name__ == '__main__':
    unittest.main()
