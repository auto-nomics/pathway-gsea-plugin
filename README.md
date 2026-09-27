# pathway-gsea plugin

Migrated from the legacy `pathway_gsea_container` wrapper in nodes-io. One
directory = one plugin family = one git-able unit. The family holds a single
node kind, `pathway_gsea`, backed by the official Bioconductor `fgsea`
package (`fgseaMultilevel`). The image provides a pinned Bioconductor R
runtime with `fgsea`, `data.table`, and `jsonlite`; gene ranks and pathway
inputs stay outside the image and are staged through the standard
`AUTONOMICS_INPUT*` / `AUTONOMICS_OUTPUT*` contract.

## Layout

- `manifest.toml` — node kind `pathway_gsea`: params, ports, resources,
  image provenance
- `scripts/pathway_gsea.R.sh` — the execution script (R source despite the
  `.sh` name, matching the family convention; `interpreter = "Rscript"` in
  the manifest governs execution), inlined by the loader at startup
- `Dockerfile` — digest-pinned `bioconductor/bioconductor:3.21-R-4.5.2`
  base + fgsea/data.table/jsonlite (moved verbatim from the repo's
  `containers/fgsea/`; that directory held no fixtures or image test
  scripts, and no Rust test referenced anything in it)
- no `[[panels]]`: the legacy wrapper bound no panels

The image runs as UID/GID 1001 (`HOME=/tmp`, `R_LIBS_USER=/tmp/R/library`)
and is compatible with isolated networking, the read-only root filesystem,
and the `/work` scratch mount.

Reference: <https://bioconductor.org/packages/fgsea/>

## Build

From this directory:

```sh
podman build --format docker -t autonomics/pathway-gsea:draft .
```

Smoke-test the installed R packages:

```sh
podman run --rm autonomics/pathway-gsea:draft \
  Rscript --vanilla -e 'cat(as.character(packageVersion("fgsea")), "\n")'
```

## Publishing

Do not publish an image as a side effect of ordinary code development.
After reviewing the build, tag the immutable manifest, push it to the GHCR
namespace, and copy the manifest digest into `image.reference` in
`manifest.toml` (it is also recorded in the workspace's
`containers/image-inventory.tsv`).

```sh
podman tag autonomics/pathway-gsea:draft \
  ghcr.io/auto-nomics/autonomics/pathway-gsea:0.1.0
podman push ghcr.io/auto-nomics/autonomics/pathway-gsea:0.1.0
podman inspect --format '{{index .RepoDigests 0}}' \
  ghcr.io/auto-nomics/autonomics/pathway-gsea:0.1.0
```

Published immutable image:

```text
ghcr.io/auto-nomics/autonomics/pathway-gsea@sha256:1eb3a32abe2f910e8a3e2c7b5cd8d9f45bf267dc7c57605103e5c30310740275
```

## Migration parity

The golden test
(`crates/container-plugin/tests/pathway_gsea_migration.rs` in the
autonomics workspace) compares the compiled `ContainerCommandSpec` against
the legacy Rust wrapper's output: image, command, outputs, resources,
timeout, and panels are equal. The script is nearly byte-identical — the
legacy wrapper already drove every parameter through `AUTONOMICS_*` env
vars, so the plugin keeps that exact env contract and only prepends the
validation guards; the fgseaMultilevel call, the GMT/TSV gene-set branch,
the report columns, and the JSON epilogue are token-for-token the legacy
embedded R code.

Deltas introduced by the plugin architecture, all following the documented
v0 patterns:

- **Kind rename**: `pathway_gsea_container` → `pathway_gsea`; the artifact
  prefix follows the kind (`/artifacts/pathway_gsea_container` →
  `/artifacts/pathway_gsea`). DAG specs referencing the old kind must be
  regenerated.
- **File-to-File ports.** The legacy node accepted optional `Any` inputs
  (DataFrames were staged to TSV by Rust code, files passed through with
  extension checks) and returned a typed DataFrame on port 0 by parsing the
  TSV report. The plugin DSL is File-to-File by policy, so both inputs are
  files (`.tsv` rank table; `.gmt` or long-format `.tsv` gene sets — the
  script always handled both) and both outputs are files (consume the
  `fgsea_report_tsv` output as a table downstream). The DataFrame
  staging/parsing convenience layer is not part of the container contract
  and is not expressible in the v0 DSL.
- **`timeout_secs` and resources are node-level constants** (3600 s;
  2 CPUs, 8Gi memory, 512 pids, 1Gi shm) instead of per-instance spec
  params. The legacy spec accepted `timeout_secs`/`cpus`/`memory` overrides
  per node instance and resolved them to these same defaults when omitted;
  the plugin DSL pins them per kind.
- **Params travel via env** (`AUTONOMICS_*`, the same variable names the
  legacy wrapper used). Numbers render as JSON literals via serde_json:
  `eps = 1e-50` renders as `"1e-50"` where the legacy Rust `f64::to_string`
  spelled out the full decimal expansion — identical after float parsing
  (the documented f64-rendering nuance). `min_size`/`max_size` render as
  plain integers, matching the legacy `usize` formatting.
- **Validation moved where the DSL can express it.** `min_size >= 1`
  (`min = 1.0`) and `eps > 0` (`exclusive_min = 0.0`) became DSL bounds.
  Non-empty column names and the cross-field `max_size >= min_size` rule
  cannot be expressed, so the script enforces them with the legacy error
  messages; failures move from node-build time to container start.
- **Slightly wider gene-set input.** The legacy wrapper required a `.gmt`
  extension for File inputs on port 1 (only staged DataFrames arrived as
  TSV); the plugin accepts a long-format `.tsv` gene-set file directly,
  which the script has always supported via its GMT/TSV branch.
