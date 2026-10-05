"""A tree block: one tree of a newick input drawn as a rectangular phylogram in inline SVG, its tips
annotated from another block. Rooted as written, ladderised, fixed row height. Every walk uses an explicit
stack or a pass over the preorder node list, so tree depth is not limited."""

import html
import math
from collections import Counter

from groundie import query as q
from groundie.query import SpecError, ident, listing, plural

ROW = 18  # px between tips
CHAR = 7.3  # px per character of 12px monospace: an estimate, there are no font metrics at build time
PAD = 8  # px around the figure
DOT = 4  # tip dot radius
GAP = 16  # px between label columns
WIDTH = 1000  # the page column: the tree takes what the labels leave of it, and at least MIN_TREE
MIN_TREE = 240
COLOURS = 3  # categorical slots in style.css; values past them share one grey


def as_dict(block):
    """The block as written, without its name."""
    return {k: v for k, v in vars(block).items() if k != "name" and v is not None}


def render(con, blocks, block):
    """The figure and legend as HTML, the lines describing it, and the manifest entry."""
    nodes, file, tree, several = read_tree(con, block)
    order, x, y, cladogram, notes = layout(nodes)
    tips = [nodes[i][2] for i in order]  # raw labels in drawing order
    repeated = sum(1 for n in Counter(t for t in tips if t is not None).values() if n > 1)
    if repeated:
        notes.append(f"{plural(repeated, 'tip label')} on more than one tip")
    where = f" ({file}, tree {tree})" if several else ""
    lines = [f"Tree from source {block.source}{where}: {plural(len(order), 'tip')}; rooted as written; "
             + "; ".join(notes)]
    entry = {"name": block.name, "block": as_dict(block), "source": block.source, "file": file, "tree": tree,
             "tips": len(order), "matched": None, "unmatched_tips": [], "unmatched_rows": [],
             "clamped_lengths": sum(1 for n in nodes[1:] if n[3] is not None and n[3] < 0)}

    rows = {}  # tip label -> the row of tips it matched, as a dict
    if block.tips:
        rows, line = match_tips(con, blocks, block, tips, entry)
        lines.append(line)
    lines.append((f"Coloured by {block.color}; labelled" if block.color else "Labelled")
                 + f" by {', '.join(block.label) if block.label else 'tip label'}")

    colour_of, legend = colours(block, tips, rows)
    svg = draw(block, nodes, order, x, y, cladogram, rows, colour_of)
    out = [svg]
    if entry["unmatched_tips"]:
        legend.append(f"grey italic{', no dot' if block.color else ''}: the tip label of a tip with no row in "
                      f"{html.escape(block.tips)}")
    if legend:
        out.append(f'<p class="legend">{"; ".join(legend)}</p>\n')
    return "".join(out), lines, entry


def read_tree(con, block):
    """The chosen tree's nodes in preorder as (node, parent, label, branch_length), its file and position,
    and whether the source holds more than one tree. The source is a newick input: checked at parse."""
    name = f"tree {block.name!r}"
    source = ident(block.source)
    trees = con.execute(f"select file, tree from {source} group by file, tree order by file, tree").fetchall()
    chosen = [t for t in trees if block.file in (None, t[0]) and block.tree in (None, t[1])]
    if len(chosen) != 1:
        asked = ", ".join(f"{k} {v}" for k, v in (("file", block.file), ("tree", block.tree)) if v is not None)
        raise SpecError(f"{name}: source {block.source!r} holds {plural(len(trees), 'tree')} "
                         f"({listing([f'{f} tree {n}' for f, n in trees])}); "
                         + (f"{len(chosen)} match {asked}; " if asked else "")
                         + "pick one with 'tree' (its position in its file, from 1) and, across files, 'file'")
    file, tree = chosen[0]
    nodes = con.execute(f"select node, parent, label, branch_length from {source} where file = ? and tree = ? "
                        "order by node", [file, tree]).fetchall()
    return nodes, file, tree, len(trees) > 1


def layout(nodes):
    """Tip indexes in drawing order, x and y per node (x in branch length units, y in rows), and notes on
    how the lengths were drawn. Nodes are in preorder: a parent comes before its children, so a forward
    pass reaches parents first and a backward pass reaches children first."""
    at = {n[0]: i for i, n in enumerate(nodes)}
    parent = [None if n[1] is None else at[n[1]] for n in nodes]
    children = [[] for _ in nodes]
    for i in range(1, len(nodes)):
        children[parent[i]].append(i)
    size = [1] * len(nodes)  # tips below each node
    for i in reversed(range(len(nodes))):
        if children[i]:
            size[i] = sum(size[c] for c in children[i])
    for kids in children:
        kids.sort(key=lambda c: size[c])  # ladderise; the sort is stable, so written order breaks ties
    order, stack = [], [0]
    while stack:
        i = stack.pop()
        if children[i]:
            stack.extend(reversed(children[i]))
        else:
            order.append(i)

    lengths = [n[3] for n in nodes[1:]]  # the root's own length is not drawn
    x = [0.0] * len(nodes)
    lengths_given = any(l is not None for l in lengths)
    if not lengths_given:
        notes = ["no branch lengths, drawn as a cladogram with tips aligned"]
        height = [0] * len(nodes)
        for i in reversed(range(len(nodes))):
            if children[i]:
                height[i] = 1 + max(height[c] for c in children[i])
        x = [float(height[0] - h) for h in height]
    else:
        notes = ["branch lengths as given"]
        for i in range(1, len(nodes)):
            x[i] = x[parent[i]] + max(nodes[i][3] or 0, 0)
        negative = sum(1 for l in lengths if l is not None and l < 0)
        missing = sum(1 for l in lengths if l is None)
        if negative:
            notes.append(f"{plural(negative, 'negative length')} drawn as zero")
        if missing:
            notes.append(f"{plural(missing, 'missing length')} drawn as zero")
    y = [0.0] * len(nodes)
    for row, i in enumerate(order):
        y[i] = float(row)
    for i in reversed(range(len(nodes))):
        if children[i]:
            y[i] = (y[children[i][0]] + y[children[i][-1]]) / 2
    return order, x, y, not lengths_given, notes


def match_tips(con, blocks, block, tips, entry):
    """The rows of the tips block keyed by match value, and the description line. Fills the
    manifest entry's match counts."""
    name = f"tree {block.name!r}"
    cur = con.execute(f"select * from {ident(block.tips)}")
    columns = [d[0] for d in cur.description]
    missing = [c for c in [block.match, *(block.label or []), block.color] if c is not None and c not in columns]
    if missing:
        raise SpecError(f"{name}: {', '.join(map(repr, missing))} not a column of {block.tips}; "
                         f"columns: {', '.join(columns)}")
    rows, shared = {}, {}
    for values in cur.fetchall():
        row = dict(zip(columns, values))
        key = None if row[block.match] is None else str(row[block.match])
        if key in rows:
            shared[key] = None
        else:
            rows[key] = row
    labels = set(tips)
    no_row = ["(no label)" if t is None else t for t in tips if t is None or t not in rows]
    no_tip = ["null" if k is None else k for k in rows if k is None or k not in labels]
    matched = len(tips) - len(no_row)
    entry.update(matched=matched, unmatched_tips=no_row, unmatched_rows=no_tip)
    line = (f"Tips from {q.about(blocks[block.tips])}; matched on {block.match}; "
            f"{matched} of {plural(len(tips), 'tip')} matched")
    if no_row:
        line += f"; tips with no row: {listing(no_row)}"
    if no_tip:
        line += f"; rows with no tip: {listing(no_tip)}"
    if shared:
        line += f"; match values on more than one row, first row used: {listing(list(shared))}"
    return {k: v for k, v in rows.items() if k is not None}, line


def colours(block, tips, rows):
    """The CSS class of each colour value, and the legend entries as HTML. The most frequent values take
    the categorical slots, first seen in drawing order breaking ties; the rest share one grey."""
    if not block.color:
        return {}, []
    values = [rows[t][block.color] for t in tips if t in rows]
    seen = Counter(str(v) for v in values if v is not None)
    first = {}
    for v in values:
        if v is not None:
            first.setdefault(str(v), len(first))
    ranked = sorted(seen, key=lambda v: (-seen[v], first[v]))
    colour_of = {v: f"c{k + 1}" if k < COLOURS else "other" for k, v in enumerate(ranked)}
    legend = [f"{html.escape(block.color)}:"]
    legend += [f"{swatch(colour_of[v])}{html.escape(v)}" for v in ranked[:COLOURS]]
    if len(ranked) > COLOURS:
        legend.append(f"{swatch('other')}other values, sharing one grey: {html.escape(listing(ranked[COLOURS:]))}")
    if any(v is None for v in values):
        legend.append(f"{swatch('null')}null")
    return colour_of, [" ".join(legend)]


def swatch(cls):
    return (f'<svg class="swatch" viewBox="0 0 10 10" width="10" height="10" aria-hidden="true">'
            f'<circle class="tree-dot {cls}" cx="5" cy="5" r="{DOT}"/></svg>')


def draw(block, nodes, order, x, y, cladogram, rows, colour_of):
    """The SVG: branches, then per tip a leader, a dot and the label columns, then the scale bar."""
    headers = block.label or []
    blank = [""] * (len(headers) or 1)

    def cells(raw):  # a tip with no row shows its raw label, set apart, in the first column
        if block.label and raw in rows:
            return ["" if rows[raw][c] is None else str(rows[raw][c]) for c in block.label]
        return ["" if raw is None else raw, *blank[1:]]

    texts = [cells(nodes[i][2]) for i in order]
    widths = [max(len(s) for s in col) * CHAR for col in zip(headers or blank, *texts)]
    labels_w = sum(widths) + GAP * (len(widths) - 1)
    tree_w = max(MIN_TREE, WIDTH - 2 * PAD - 3 * DOT - 12 - labels_w)
    depth = max(x)
    per_unit = tree_w / depth if depth > 0 else 0
    label_x = PAD + tree_w + 3 * DOT + 12
    top = PAD + ROW * (1.5 if headers else 0.5)
    bar = 0 if cladogram or depth == 0 else nice(depth * 0.3)
    width = label_x + labels_w + PAD
    height = top + (len(order) - 1) * ROW + ROW / 2 + (28 if bar else PAD)

    def px(i):
        return PAD + x[i] * per_unit

    def py(i):
        return top + y[i] * ROW

    out = [f'<svg class="tree" viewBox="0 0 {width:.0f} {height:.0f}" width="{width:.0f}" height="{height:.0f}" '
           f'role="img" aria-label="Tree {html.escape(block.name)}">\n']
    at = {n[0]: i for i, n in enumerate(nodes)}
    for i, n in enumerate(nodes[1:], 1):
        p = at[n[1]]
        out.append(f'<path class="branch" d="M{px(p):.1f} {py(p):.1f}V{py(i):.1f}H{px(i):.1f}"/>\n')
    col_x = [label_x + sum(widths[:k]) + GAP * k for k in range(len(widths))]
    for k, h in enumerate(headers):
        out.append(f'<text class="head" x="{col_x[k]:.1f}" y="{top - ROW + 4:.1f}">{html.escape(h)}</text>\n')
    for i, row_texts in zip(order, texts):
        raw = nodes[i][2]
        cy = py(i)
        out.append(f'<g class="tip"><title>{html.escape("(no label)" if raw is None else raw)}</title>')
        start = px(i) + (DOT + 3 if block.color else 3)
        if label_x - 4 > start:
            out.append(f'<path class="leader" d="M{start:.1f} {cy:.1f}H{label_x - 4:.1f}"/>')
        row = rows.get(raw)
        if block.color and row is not None:
            v = row[block.color]
            cls = "null" if v is None else colour_of[str(v)]
            out.append(f'<circle class="tree-dot {cls}" cx="{px(i):.1f}" cy="{cy:.1f}" r="{DOT}"/>')
        muted = ' class="unmatched"' if block.tips and row is None else ""
        for k, s in enumerate(row_texts):
            if s:
                out.append(f'<text{muted} x="{col_x[k]:.1f}" y="{cy + 4:.1f}">{html.escape(s)}</text>')
        out.append("</g>\n")
    if bar:
        by = height - 14
        end = PAD + bar * per_unit
        out.append(f'<path class="scale" d="M{PAD} {by - 3:.1f}v6M{PAD} {by:.1f}H{end:.1f}M{end:.1f} {by - 3:.1f}v6"/>\n'
                   f'<text class="scale" x="{end + 6:.1f}" y="{by + 4:.1f}">{bar:g}</text>\n')
    out.append("</svg>\n")
    return "".join(out)


def nice(target):
    """The largest of 1, 2 or 5 times a power of ten that is at most target."""
    step = 10.0 ** math.floor(math.log10(target))
    return max(m * step for m in (1, 2, 5) if m * step <= target)
