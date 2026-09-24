# Configuration

pyLocusZoom is a Python library, so plot options are keyword arguments and
frozen config models passed to each plot call. Those are documented once, in
the [User Guide](USER_GUIDE.md#plotter-reference):
[`ColumnConfig`](USER_GUIDE.md#columnconfig),
[`DisplayConfig`](USER_GUIDE.md#displayconfig),
[`LDConfig`](USER_GUIDE.md#ldconfig),
[`LiftoverConfig`](USER_GUIDE.md#liftoverconfig),
[`PanelInputs`](USER_GUIDE.md#panelinputs),
[`GenomeWideConfig`](USER_GUIDE.md#genomewideconfig),
[`GenomeWideStyle`](USER_GUIDE.md#genomewidestyle) and
[`ColocConfig`](USER_GUIDE.md#colocconfig).

This page covers what is configured outside a plot call: the on-disk caches for
reference data, and the environment variables that move them.

There is no `.env` file and no runtime configuration file. The library does
not read any `PYLOCUSZOOM_*` environment variables.

## Cache Location

Recombination maps, liftover chains and gene annotations download once and are
cached under one base directory:

| Platform | Cache base directory |
| -------- | -------------------- |
| Linux, macOS | `$XDG_CACHE_HOME/pylocuszoom`, or `~/.cache/pylocuszoom` when `XDG_CACHE_HOME` is unset |
| Windows | `%LOCALAPPDATA%\pylocuszoom`, or `~\AppData\Local\pylocuszoom` when `LOCALAPPDATA` is unset |
| Databricks (any machine where `/dbfs` exists) | `/dbfs/FileStore/reference_data`, whatever `XDG_CACHE_HOME` says |

Each kind of data has its own subfolder of the base:

| Subfolder | Contents |
| --------- | -------- |
| `recombination_maps/` | The managed recombination maps, downloaded on the first plot that needs them (canine only) |
| `liftover/` | Liftover chain files, such as the CanFam3.1 to CanFam4 chain the recombination maps are lifted through |
| `ensembl/{ensembl_species}/` | Genes and exons fetched from Ensembl, one ZIP per region; the folder is Ensembl's species name, such as `homo_sapiens` |
| `ucsc/{ucsc_genome}/` | Genes and exons fetched from UCSC for CanFam3.1, CanFam4 and FelCat9, such as `ucsc/canFam3/` |

Replacing the map set never touches the chains, and a region fetched from
Ensembl never collides with the same region fetched from UCSC.
`clear_gene_cache("ensembl")` or `clear_gene_cache("ucsc")` empties one gene
cache; see [Automatic Gene Annotations](USER_GUIDE.md#automatic-gene-annotations-from-ensembl).

The base directory is resolved in one place,
[`utils._platform_cache_base()`](../src/pylocuszoom/utils.py).
[`recombination.get_default_data_dir()`](../src/pylocuszoom/recombination.py)
returns the map folder, and
[`_gene_cache.cache_root()`](../src/pylocuszoom/_gene_cache.py) the gene folders.

## Environment Variables

These are the only environment variables the library reads:

| Variable | Platform | Effect |
| -------- | -------- | ------ |
| `XDG_CACHE_HOME` | Linux, macOS | Moves the cache base to `$XDG_CACHE_HOME/pylocuszoom`. Ignored on Windows and on Databricks. |
| `LOCALAPPDATA` | Windows | Moves the cache base to `%LOCALAPPDATA%\pylocuszoom`. |

## Your Own Recombination Maps

To pre-download maps into a chosen directory, call
`download_canine_recombination_maps(output_dir="/path/to/maps")`. The directory
must be new, empty or hold only a previous map set; one holding other files
raises `ValidationError` and is left untouched.

Passing `recomb_data_dir` to the plotter, or `data_dir` to the map helpers,
selects a read-only caller directory. It never downloads or replaces files there,
and its coordinates must already use the requested build. With no directory,
`ensure_recomb_maps()` manages the platform cache and may download built-in maps.
The file format is in
[Recombination Map Files](USER_GUIDE.md#recombination-map-files).

## Per-Environment Overrides

pyLocusZoom does not distinguish "development" vs "production" environments
at runtime — it is a library, not a service. There are no
`.env.development` / `.env.production` files and no `NODE_ENV`-style switch.

If you need per-environment behaviour, do it at the caller level, e.g.:

- Set `XDG_CACHE_HOME` / `LOCALAPPDATA` per machine to control where
  reference data is cached.
- Pre-download canine maps with `download_canine_recombination_maps(output_dir=...)`,
  then pass that directory as `recomb_data_dir` to the plotter. Caller maps must
  already use the requested genome build.
- On Databricks, the `/dbfs/FileStore/reference_data` base is selected
  automatically.

The settings in `pyproject.toml` (pytest, ruff, coverage) affect only
contributors; [DEVELOPMENT.md](DEVELOPMENT.md) covers them.
