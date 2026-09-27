# Official Bioconductor fgsea run driven by the plugin env channel: every
# legacy wrapper spec param arrives as an environment variable instead of
# being interpolated into R source by Rust. (The legacy wrapper already
# used this exact env contract, so the body below is byte-identical to the
# embedded legacy R script.) Numbers render as JSON literals;
# as.integer()/as.numeric() parse them exactly like the legacy embedded
# literals — eps renders as "1e-50", the same double the legacy Rust
# decimal spelling produced.

gene_col <- Sys.getenv("AUTONOMICS_GENE_COL")
score_col <- Sys.getenv("AUTONOMICS_SCORE_COL")
pathway_col <- Sys.getenv("AUTONOMICS_PATHWAY_COL")
pathway_gene_col <- Sys.getenv("AUTONOMICS_PATHWAY_GENE_COL")
min_size <- as.integer(Sys.getenv("AUTONOMICS_MIN_SIZE"))
max_size <- as.integer(Sys.getenv("AUTONOMICS_MAX_SIZE"))
eps <- as.numeric(Sys.getenv("AUTONOMICS_EPS"))

# The manifest DSL cannot express non-empty strings or cross-field
# comparisons (the legacy Rust validate()); the script enforces them with
# the legacy error messages. min_size >= 1 and eps > 0 are covered by the
# DSL bounds (min = 1.0, exclusive_min = 0.0).
if (!nzchar(trimws(gene_col))) stop("gene_col cannot be empty")
if (!nzchar(trimws(score_col))) stop("score_col cannot be empty")
if (!nzchar(trimws(pathway_col))) stop("pathway_col cannot be empty")
if (!nzchar(trimws(pathway_gene_col))) stop("pathway_gene_col cannot be empty")
if (max_size < min_size) stop("max_size must be at least min_size")

rank_table <- data.table::fread(Sys.getenv("AUTONOMICS_INPUT0"), check.names = FALSE)
if (!gene_col %in% names(rank_table)) stop("rank table has no column: ", gene_col)
if (!score_col %in% names(rank_table)) stop("rank table has no column: ", score_col)
genes <- as.character(rank_table[[gene_col]])
if (anyNA(genes) || any(!nzchar(genes))) stop("rank gene column contains null or empty values")
if (anyDuplicated(genes)) stop("rank gene column contains duplicate values")
scores <- as.numeric(rank_table[[score_col]])
if (anyNA(scores)) stop("rank score column contains null or nonnumeric values")
stats <- setNames(scores, genes)

if (grepl("\\.gmt$", Sys.getenv("AUTONOMICS_INPUT1"))) {
  gmt_lines <- readLines(Sys.getenv("AUTONOMICS_INPUT1"), warn = FALSE)
  gmt_lines <- gmt_lines[nzchar(gmt_lines)]
  if (!length(gmt_lines)) stop("GMT file contains no gene sets")
  fields <- strsplit(gmt_lines, "\t", fixed = TRUE)
  if (any(lengths(fields) < 3L)) stop("GMT contains a set without genes")
  pathways <- lapply(fields, function(field) field[-c(1, 2)])
  names(pathways) <- vapply(fields, `[`, character(1), 1L)
  if (anyDuplicated(names(pathways))) stop("GMT contains duplicate pathway names")
} else {
  pathway_table <- data.table::fread(Sys.getenv("AUTONOMICS_INPUT1"), check.names = FALSE)
  if (!pathway_col %in% names(pathway_table)) stop("pathway table has no column: ", pathway_col)
  if (!pathway_gene_col %in% names(pathway_table)) stop("pathway table has no column: ", pathway_gene_col)
  set_names <- as.character(pathway_table[[pathway_col]])
  set_genes <- as.character(pathway_table[[pathway_gene_col]])
  if (anyNA(set_names) || any(!nzchar(set_names))) stop("pathway column contains null or empty values")
  if (anyNA(set_genes) || any(!nzchar(set_genes))) stop("pathway gene column contains null or empty values")
  pathways <- split(set_genes, set_names)
}
pathways <- lapply(pathways, unique)

result <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = min_size,
  maxSize = max_size,
  eps = eps
)
result$leadingEdge <- vapply(result$leadingEdge, paste, character(1), collapse = ";")
report <- result[, .(pathway, pval, padj, NES, size, leadingEdge)]
data.table::fwrite(report, Sys.getenv("AUTONOMICS_OUTPUT0"), sep = "\t", quote = FALSE)

json <- list(
  method = "fgseaMultilevel",
  parameters = list(
    min_size = min_size,
    max_size = max_size,
    eps = eps,
    gene_col = gene_col,
    score_col = score_col,
    pathway_col = pathway_col,
    pathway_gene_col = pathway_gene_col
  ),
  results = report
)
writeLines(
  jsonlite::toJSON(json, auto_unbox = TRUE, pretty = TRUE, na = "null"),
  Sys.getenv("AUTONOMICS_OUTPUT1")
)
