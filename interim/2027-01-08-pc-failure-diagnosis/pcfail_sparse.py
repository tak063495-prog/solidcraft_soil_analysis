"""Resolve the existing CI sparse OR conditions without OS-dependent paths."""
import re


def rules(workflow,group):
    text=workflow.read_text(encoding='utf-8')
    block=text.split('          sparse-checkout: |\n',1)[1].split('      - name:',1)[0]
    output=[]
    for raw in block.splitlines():
        line=raw.strip()
        if not line:continue
        if line.startswith('${{'):
            match=re.fullmatch(r"\$\{\{ (.*?) && '([^']*)' \|\| '([^']*)' \}\}",line)
            if not match:raise ValueError('unreviewed CI conditional')
            condition,yes,no=match.groups();names=re.findall(r"matrix.group == '([^']+)'",condition)
            remainder=re.sub(r"matrix.group == '[^']+'",'',condition)
            if not names or not re.fullmatch(r'[()\s|]*',remainder):raise ValueError('unreviewed CI group expression')
            output.append(yes if group in names else no)
        else:output.append(line)
    return output
