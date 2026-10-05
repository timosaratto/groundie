"""Newick trees, one row per node. Several trees per file and several files per source are allowed.
Facts are kept raw: a label is unquoted but otherwise as written (underscores stay), a missing length is
null, a negative length stays negative, comments are kept as written, brackets included. A numeric label
on an internal node is also stored as 'support'. Nodes are numbered in preorder from 1 within each tree.
The parser is iterative, so tree depth is not limited by Python's recursion limit."""

import re
from pathlib import Path

NAME = "newick"
VERSION = "0.1"
OPTIONS = ()
COLUMNS = (("file", "varchar"), ("tree", "integer"), ("node", "integer"), ("parent", "integer"),
           ("label", "varchar"), ("branch_length", "double"), ("support", "double"),
           ("is_tip", "boolean"), ("comment", "varchar"))

NUMBER = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
DELIMITERS = set("()[]',:;") | set(" \t\r\n")


def rows(path, options):
    path = Path(path)
    return [(path.name, *node) for node in parse(path.read_text(), path.name)]


def tokens(text, name):
    """Yield (position, kind, value): kind is one of ( ) , : ; label comment."""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            i += 1
        elif c in "(),:;":
            yield i, c, c
            i += 1
        elif c == "[":
            end = text.find("]", i)
            if end < 0:
                raise ValueError(f"{name}, character {i + 1}: comment '[' is never closed")
            yield i, "comment", text[i : end + 1]
            i = end + 1
        elif c == "]":
            raise ValueError(f"{name}, character {i + 1}: ']' without an opening '['")
        elif c == "'":
            start, i, parts = i, i + 1, []
            while True:
                end = text.find("'", i)
                if end < 0:
                    raise ValueError(f"{name}, character {start + 1}: quoted label is never closed")
                parts.append(text[i:end])
                if text.startswith("''", end):
                    parts.append("'")
                    i = end + 2
                else:
                    i = end + 1
                    break
            yield start, "label", "".join(parts)
        else:
            start = i
            while i < n and text[i] not in DELIMITERS:
                i += 1
            yield start, "label", text[start:i]


def parse(text, name):
    """Rows of (tree, node, parent, label, branch_length, support, is_tip, comment) for every tree in text."""
    out = []  # rows as lists, finished into tuples at the end
    tree, stack, current = 0, [], None  # stack: open internal nodes (row indexes); current: the node being written
    expecting = True  # at the start of a tree, after '(' and after ',': the next token begins a node
    started = False  # the current tree has content
    pending_comment = ""  # comment seen where a node is expected: it belongs to that node
    after_colon = False

    first = 0  # row index of the current tree's first node

    def new_node(is_tip):
        nonlocal tree, first, current, started, pending_comment
        if not started:
            tree, first, started = tree + 1, len(out), True
        parent = out[stack[-1]][1] if stack else None
        out.append([tree, len(out) - first + 1, parent, None, None, None, is_tip, pending_comment or None])
        pending_comment = ""
        current = len(out) - 1

    for pos, kind, value in tokens(text, name):
        if after_colon:
            if kind != "label" or not NUMBER.fullmatch(value):
                raise ValueError(f"{name}, character {pos + 1}: expected a branch length after ':'")
            out[current][4] = float(value)
            after_colon = False
            continue
        if kind == "comment":
            if expecting:
                pending_comment += value
            else:
                out[current][7] = (out[current][7] or "") + value
            continue
        if kind == "(":
            if not expecting:
                raise ValueError(f"{name}, character {pos + 1}: unexpected '('")
            new_node(False)
            stack.append(current)
            continue
        if expecting and kind in ",):":
            new_node(True)  # a tip with no label
        if expecting and kind == "label":
            new_node(True)
            out[current][3] = value
            expecting = False
        elif kind == "label":
            row = out[current]
            if row[3] is not None or row[4] is not None:
                raise ValueError(f"{name}, character {pos + 1}: unexpected label {value!r}")
            row[3] = value
        elif kind == ":":
            if out[current][4] is not None:
                raise ValueError(f"{name}, character {pos + 1}: a second branch length")
            expecting, after_colon = False, True
        elif kind == ",":
            if not stack:
                raise ValueError(f"{name}, character {pos + 1}: ',' outside parentheses")
            expecting = True
        elif kind == ")":
            if not stack:
                raise ValueError(f"{name}, character {pos + 1}: unbalanced ')'")
            current = stack.pop()
            expecting = False
        elif kind == ";":
            if stack:
                raise ValueError(f"{name}, character {pos + 1}: unbalanced '(': {len(stack)} not closed")
            if not started:
                raise ValueError(f"{name}, character {pos + 1}: empty tree before ';'")
            started, expecting, current = False, True, None
    if stack:
        raise ValueError(f"{name}, character {len(text)}: unbalanced '(': {len(stack)} not closed")
    if after_colon:
        raise ValueError(f"{name}, character {len(text)}: expected a branch length after ':'")
    if started or pending_comment:
        raise ValueError(f"{name}, character {len(text)}: missing ';' at the end of the last tree")
    if not out:
        raise ValueError(f"{name}: no tree found")
    for row in out:
        if not row[6] and row[3] is not None and NUMBER.fullmatch(row[3]):
            row[5] = float(row[3])
    return [tuple(row) for row in out]
