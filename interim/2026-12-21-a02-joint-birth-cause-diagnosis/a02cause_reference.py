"""Independent diagnostic partition policy, isolated from the frozen reference.

Only the initial contact-partition count changes; Newton, forces, fields,
material checks and maximum-increment enforcement remain the frozen algorithm.
The default minimum of one retains the original numerical operations.
"""
import ast
import inspect
from a02cause_binding import reference


def controlled_run(minimum_parts=1):
    if type(minimum_parts) is not int or not 1 <= minimum_parts <= 16:
        raise ValueError('diagnostic minimum partitions must be integer 1..16')
    if minimum_parts == 1:return reference.run
    tree=ast.parse(inspect.getsource(reference.run));changed=0
    for node in ast.walk(tree):
        if (isinstance(node,ast.Assign) and len(node.targets)==1
            and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='counts'
            and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute)
            and isinstance(node.value.func.value,ast.Name) and node.value.func.value.id=='np'
            and node.value.func.attr=='ones'):
            assert len(node.value.args)==1 and isinstance(node.value.args[0],ast.Name) and node.value.args[0].id=='nq'
            node.value.func.attr='full';node.value.args.append(ast.Constant(minimum_parts));changed+=1
    if changed != 1:raise ValueError('frozen reference partition initialization changed')
    ast.fix_missing_locations(tree);namespace=dict(vars(reference))
    exec(compile(tree,'a02cause_minimum_parts_'+str(minimum_parts),'exec'),namespace)
    return namespace['run']
