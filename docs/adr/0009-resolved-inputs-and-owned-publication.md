# ADR 0009: Resolve data roles before rendering and publish only owned state

Status: Accepted

## Context

The figure and panel model already gives every plot one rendering path. Its
inputs still allowed several interpretations downstream: a numeric lead position
could identify a different row, merged column suffixes could shadow configured
fields, and a temporary alias resolution could be lost by a composed plot.
Filesystem paths similarly failed to distinguish caller data from managed caches
or an incomplete write from a published entry.

## Decision

Resolve each regional frame's selected rows, column names, LD options and lead
identity before panel construction. Genome-wide preparation projects configured
roles per frame before building shared layout. Colocalization projects each
source's requested roles before merging; internal column names have one meaning.
Chromosome and absolute position define coordinate overlap. Position-only callers
must explicitly identify their shared chromosome. Fine-mapping loaders preserve
reported PIPs and membership without making an implicit statistical inference.

Use the existing canonical heatmap edges in drawing adapters and SNP outlines.
Order regional heatmap coordinates with both matrix axes before drawing them.
Keep FigurePlan, panel-owned drawing, the shared layout and public config values.

An explicit map directory is caller-owned and read-only. Only the default managed
cache may install a built-in map set or apply its known assembly conversion.
Map archives supply regular text members, streamed under generated canonical
names. General filesystem extraction is unnecessary.

Each HTTP download owns a private temporary file. Gene and exon annotations form
one cache entry, an archive replaced only after both members have been written.
A failed update leaves the previous entry available; unreadable entries are
cache misses. Example verification likewise generates outside the checkout and
only copies outputs back through an explicit acceptance operation.

## Consequences

Extra input metadata cannot redefine selected column roles. One prepared regional
input supplies all association consumers. Writers cannot share an in-progress
filename or publish half an annotation entry. Custom map files survive plotting.

Position-only overlap calls must supply `common_chrom`. FINEMAP/CAVIAR calls that
relied on inferred membership now receive PIPs without a credible-set column.
Legacy annotation CSV pairs are cold misses. These behavior changes are deliberate
corrections to ambiguous or scientifically incorrect contracts; migration details
are in the user guide and changelog.

Regression tests observe native plotted points and cell bounds, loaded data,
subprocess paths resolved from its working directory, concurrent publication and
preservation of caller files. The example comparison checks serialized exports
and returns nonzero on a difference.
