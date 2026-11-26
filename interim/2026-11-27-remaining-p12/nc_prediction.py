"""Enclose an unresolved p10 sign bracket without assigning it a successful check."""
import inspect
import json


from nc_paths import predict as original


def build_predictor():
    # Keep the p10 grid, cuts, iteration count, point checks and width expression verbatim.
    # A failed midpoint leaves a *closed interval* for every possible p10 left endpoint.
    src=inspect.getsource(original)
    def insert(old,new):
        nonlocal src
        assert src.count(old)==1
        src=src.replace(old,new)
    insert("    roots = []","    roots = []\n    unresolved = []")
    insert('                mid=sample((a["width"]+b["width"])/2)', '''                midpoint=(a["width"]+b["width"])/2
                try:
                    mid=sample(midpoint)
                except AssertionError as exc:
                    failure=dict(point=i,left=a,right=b,failed_midpoint=midpoint,diagnostic=str(exc))
                    unresolved.append(failure)
                    print(json.dumps(dict(parent=index,load=lam,unresolved=failure)),file=local)
                    break''')
    insert("            roots.append(", "            else:\n                roots.append(")
    insert("    return row",'''    row["unresolved_brackets"]=unresolved
    # A partial bracket may make the width indeterminate. Never advance that parent.
    row["unique_width"]=all(.9*r["left"]["width"]>=proposed for r in unresolved)
    if not row["unique_width"]:
        row["upper_candidate_m"]=proposed
        row["predicted"]=None
    return row''')
    ns=original.__globals__.copy()
    exec(compile(src,"p12-p10-width-prediction","exec"),ns)
    return ns["predict"]


predict=build_predictor()
