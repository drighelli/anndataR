#' Convert an AnnData object to a SpatialExperiment
#'
#' @param adata The AnnData object to convert
#' @param x_mapping A string specifying the assay name for `adata$X`
#' @param assays_mapping See [as_SingleCellExperiment()]
#' @param colData_mapping See [as_SingleCellExperiment()]
#' @param rowData_mapping See [as_SingleCellExperiment()]
#' @param reducedDims_mapping See [as_SingleCellExperiment()]
#' @param colPairs_mapping See [as_SingleCellExperiment()]
#' @param rowPairs_mapping See [as_SingleCellExperiment()]
#' @param metadata_mapping See [as_SingleCellExperiment()]
#' @param spatialCoords_mapping Name of the key in `adata$obsm` to import as
#' spatial coordinates. Set to `FALSE` to skip.
#' @param imgData_mapping Logical; whether to import `adata$uns[["spatial"]]`
#' as `imgData`.
#' @param library_id One of `"auto"`, `"sample_id"`, or a column name in
#' `colData` to use as sample/library identifier.
#'
#' @return A SpatialExperiment object
#' @noRd
as_SpatialExperiment <- function(
    adata,
    x_mapping=NULL,
    assays_mapping=TRUE,
    colData_mapping=TRUE,
    rowData_mapping=TRUE,
    reducedDims_mapping=TRUE,
    colPairs_mapping=TRUE,
    rowPairs_mapping=TRUE,
    metadata_mapping=TRUE,
    spatialCoords_mapping="spatial",
    imgData_mapping=TRUE,
    library_id="auto"
) {
    check_requires(
        "Converting AnnData to SpatialExperiment",
        "SpatialExperiment",
        "Bioc"
    )

    if (!(inherits(adata, "AbstractAnnData"))) {
        cli_abort(
            "{.arg adata} must be a {.cls AbstractAnnData} but has class {.cls {class(adata)}}"
        )
    }

    if (!identical(spatialCoords_mapping, FALSE) && !is.character(spatialCoords_mapping)) {
        cli_abort(
            "{.arg spatialCoords_mapping} must be a single character string or FALSE"
        )
    }

    if (isTRUE(reducedDims_mapping) && !identical(spatialCoords_mapping, FALSE)) {
        reducedDims_mapping <- .as_SPE_drop_spatial_from_reducedDims_guess(
            adata=adata,
            spatial_key=spatialCoords_mapping[[1]]
        )
    }

    if (isTRUE(metadata_mapping)) {
        metadata_mapping <- .as_SPE_drop_reserved_uns_guess(adata)
    }

    sce <- as_SingleCellExperiment(
        adata=adata,
        x_mapping=x_mapping,
        assays_mapping=assays_mapping,
        colData_mapping=colData_mapping,
        rowData_mapping=rowData_mapping,
        reducedDims_mapping=reducedDims_mapping,
        colPairs_mapping=colPairs_mapping,
        rowPairs_mapping=rowPairs_mapping,
        metadata_mapping=metadata_mapping
    )

    coords <- NULL
    if (!identical(spatialCoords_mapping, FALSE)) {
        coords <- .as_SPE_extract_spatialCoords(
            adata=adata,
            spatial_key=spatialCoords_mapping[[1]]
        )
    }

    spe <- SpatialExperiment::SpatialExperiment(
        assays=SummarizedExperiment::assays(sce),
        rowData=SummarizedExperiment::rowData(sce),
        colData=SummarizedExperiment::colData(sce),
        metadata=S4Vectors::metadata(sce),
        spatialCoords=coords
    )

    if (length(SingleCellExperiment::reducedDims(sce)) > 0L) {
        SingleCellExperiment::reducedDims(spe) <- SingleCellExperiment::reducedDims(sce)
    }

    for (nm in SingleCellExperiment::colPairNames(sce)) {
        SingleCellExperiment::colPair(spe, nm) <- SingleCellExperiment::colPair(sce, nm)
    }

    for (nm in SingleCellExperiment::rowPairNames(sce)) {
        SingleCellExperiment::rowPair(spe, nm) <- SingleCellExperiment::rowPair(sce, nm)
    }

    if (isTRUE(imgData_mapping)) {
        img_data <- .as_SPE_extract_imgData(
            adata=adata,
            col_data=SummarizedExperiment::colData(spe),
            library_id=library_id
        )

        if (!is.null(img_data) && nrow(img_data) > 0L) {
            SpatialExperiment::imgData(spe) <- img_data
        }
    }

    spe
}


.as_SPE_drop_spatial_from_reducedDims_guess <- function(adata, spatial_key="spatial") {
    mapping <- .as_SCE_guess_reducedDims(adata)

    if (length(mapping) == 0L) {
        return(mapping)
    }

    mapping[names(mapping) != spatial_key]
}


.as_SPE_drop_reserved_uns_guess <- function(adata) {
    mapping <- .as_SCE_guess_all(adata, slot="uns")

    if (length(mapping) == 0L) {
        return(mapping)
    }

    reserved <- c("spatial", "anndataR_spatialexperiment")
    mapping[!(mapping %in% reserved)]
}


.as_SPE_get_roundtrip_uns <- function(adata) {
    if ("anndataR_spatialexperiment" %in% adata$uns_keys()) {
        return(adata$uns[["anndataR_spatialexperiment"]])
    }
    NULL
}


.as_SPE_extract_spatialCoords <- function(adata, spatial_key="spatial") {
    if (!(spatial_key %in% adata$obsm_keys())) {
        return(NULL)
    }

    coords <- adata$obsm[[spatial_key]]
    coords <- as.matrix(coords)

    if (!is.numeric(coords)) {
        cli_abort(
            "{.code adata$obsm[[{spatial_key}]]} must be numeric to be imported as spatial coordinates"
        )
    }

    roundtrip <- .as_SPE_get_roundtrip_uns(adata)
    coord_names <- NULL
    if (!is.null(roundtrip) && "spatialCoords_colnames" %in% names(roundtrip)) {
        coord_names <- unlist(roundtrip[["spatialCoords_colnames"]], use.names=FALSE)
    }

    if (is.null(coord_names) || length(coord_names) != ncol(coords)) {
        if (ncol(coords) == 2L) {
            coord_names <- c("x", "y")
        } else {
            coord_names <- paste0("spatial_", seq_len(ncol(coords)))
        }
    }

    colnames(coords) <- coord_names
    rownames(coords) <- adata$obs_names

    coords
}


.as_SPE_guess_obs_library_ids <- function(col_data, library_id="auto") {
    n <- nrow(col_data)

    if (identical(library_id, "auto")) {
        if ("sample_id" %in% names(col_data)) {
            return(as.character(col_data[["sample_id"]]))
        }
        return(rep("library1", n))
    }

    if (identical(library_id, "sample_id")) {
        if (!("sample_id" %in% names(col_data))) {
            cli_abort("Column {.val sample_id} not found in imported {.code colData}")
        }
        return(as.character(col_data[["sample_id"]]))
    }

    if (is.character(library_id) && length(library_id) == 1L) {
        if (!(library_id %in% names(col_data))) {
            cli_abort(
                "Column {.val {library_id}} not found in imported {.code colData}"
            )
        }
        return(as.character(col_data[[library_id]]))
    }

    cli_abort(
        "{.arg library_id} must be one of {.val auto}, {.val sample_id}, or a valid column name"
    )
}


.as_SPE_extract_imgData <- function(adata, col_data, library_id="auto") {
    if (!("spatial" %in% adata$uns_keys())) {
        return(NULL)
    }

    spatial_uns <- adata$uns[["spatial"]]

    if (is.null(spatial_uns) || length(spatial_uns) == 0L) {
        return(NULL)
    }

    out <- list()
    idx <- 1L

    for (lib in names(spatial_uns)) {
        x <- spatial_uns[[lib]]

        images <- if ("images" %in% names(x)) x[["images"]] else list()
        scalefactors <- if ("scalefactors" %in% names(x)) x[["scalefactors"]] else list()
        metadata <- if ("metadata" %in% names(x)) x[["metadata"]] else list()

        image_names <- union(names(images), names(metadata))
        if (length(image_names) == 0L) {
            image_names <- "image_1"
        }

        for (img_nm in image_names) {
            row_i <- list(
                sample_id=lib,
                image_id=img_nm
            )

            if (img_nm %in% names(images)) {
                row_i[["data"]] <- list(images[[img_nm]])
            } else {
                row_i[["data"]] <- list(NULL)
            }

            if (length(scalefactors) > 0L) {
                for (nm in names(scalefactors)) {
                    row_i[[nm]] <- scalefactors[[nm]]
                }
            }

            if (img_nm %in% names(metadata)) {
                md_i <- metadata[[img_nm]]
                if (is.list(md_i) && length(md_i) > 0L) {
                    for (nm in names(md_i)) {
                        row_i[[nm]] <- md_i[[nm]]
                    }
                }
            }

            out[[idx]] <- row_i
            idx <- idx + 1L
        }
    }

    if (length(out) == 0L) {
        return(NULL)
    }

    out_df <- S4Vectors::DataFrame(do.call(rbind, lapply(out, function(x) {
        as.data.frame(x, stringsAsFactors=FALSE)
    })))

    list_cols <- vapply(out, function(x) "data" %in% names(x), logical(1))
    if (any(list_cols)) {
        out_df[["data"]] <- S4Vectors::I(lapply(out, `[[`, "data"))
        out_df[["data"]] <- lapply(out_df[["data"]], function(x) x[[1]])
    }

    out_df
}
