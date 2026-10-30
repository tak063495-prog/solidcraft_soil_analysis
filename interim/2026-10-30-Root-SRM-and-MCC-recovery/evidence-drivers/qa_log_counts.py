"""Strict unittest inventory; buffered application stdout may follow the summary."""
import re

def unittest_summary(text, expected):
    summaries = re.findall(r'^Ran (\d+) tests in ([\d.]+)s$', text, re.MULTILINE)
    assert len(summaries) == 1 and int(summaries[0][0]) == expected
    assert len(re.findall(r'^OK$', text, re.MULTILINE)) == 1
    starts = list(re.finditer(r'^test_\S+ \(([^)\n]+)\) \.\.\. ', text, re.MULTILINE))
    cases = [m.group(1) for m in starts]
    assert len(cases) == len(set(cases)) == expected
    summary_start = re.search(r'^-{5,}\nRan \d+ tests in ', text, re.MULTILINE)
    assert summary_start is not None
    interleaved = 0
    for index, start in enumerate(starts):
        stop = starts[index + 1].start() if index + 1 < len(starts) else summary_start.start()
        block = text[start.start():stop]
        lines = [line.rstrip() for line in block.splitlines() if line.strip()]
        assert lines
        if lines[-1] == 'ok':
            interleaved += 1
        else:
            assert re.fullmatch(r'test_\S+ \([^)\n]+\) \.\.\. ok', lines[-1])
    assert not re.search(r'^FAILED\b|\.\.\. (?:FAIL|ERROR|skipped|expected failure|unexpected success)\b', text, re.MULTILINE)
    return dict(passed=True, tests=expected, unique_successful_cases=len(cases), failed=0, skipped=0,
        actual_elapsed_seconds=float(summaries[0][1]),
        single_exact_summary_and_OK_and_all_unique_test_cases_verified=True,
        test_successes_after_interleaved_tool_output=interleaved,
        trailing_buffered_application_stdout_not_treated_as_a_test_result=True)
