#' Convert a SpatialExperiment object to an AnnData object
#'
#' @param spe A SpatialExperiment::SpatialExperiment object
#' @param x_mapping See [as_AnnData()]
#' @param layers_mapping See [as_AnnData()]
#' @param obs_mapping See [as_AnnData()]
#' @param var_mapping See [as_AnnData()]
#' @param obsm_mapping See [as_AnnData()]
#' @param varm_mapping See [as_AnnData()]
#' @param obsp_mapping See [as_AnnData()]
#' @param varp_mapping See [as_AnnData()]
#' @param uns_mapping See [as_AnnData()]
#' @param spatialCoords_mapping Name of the key in `adata$obsm` where
#' spatial coordinates will be stored. Set to `FALSE` to skip exporting
#' spatial coordinates.
#' @param imgData_mapping Logical; whether to export `imgData(spe)` to
#' `adata$uns[["spatial"]]`.
#' @param library_id One of `"auto"`, `"sample_id"`, `"sample_idImg"`, or a
#' single character string naming a column in `colData(spe)` to use as the
#' library identifier.
#' @param image_mode One of `"embed"`, `"reference"`, or `"none"`.
#' @param output_class See [as_AnnData()]
#' @param ... See [as_AnnData()]
#'
#' @return An AnnData object
#' @export
#' @importFrom SpatialExperiment SpatialExperiment spatialCoords imgData
#' @noRd
from_SpatialExperiment <- function(
    spe,
    x_mapping=TRUE,
    layers_mapping=TRUE,
    obs_mapping=TRUE,
    var_mapping=TRUE,
    obsm_mapping=TRUE,
    varm_mapping=TRUE,
    obsp_mapping=TRUE,
    varp_mapping=TRUE,
    uns_mapping=TRUE,
    spatialCoords_mapping="spatial",
    imgData_mapping=TRUE,
    library_id="auto",
    image_mode=c("embed", "reference", "none"),
    output_class=c("InMemory", "HDF5AnnData", "ReticulateAnnData"),
    ...
) {
    check_requires(
        "Converting SpatialExperiment to AnnData",
        "SpatialExperiment",
        "Bioc"
    )

    image_mode <- match.arg(image_mode)
    output_class <- match.arg(output_class)

    if (!(inherits(spe, "SpatialExperiment"))) {
        cli_abort(
            "{.arg spe} must be a {.cls SpatialExperiment} but has class {.cls {class(spe)}}"
        )
    }

    if (!identical(spatialCoords_mapping, FALSE) && !is.character(spatialCoords_mapping)) {
        cli_abort(
            "{.arg spatialCoords_mapping} must be a single character string or FALSE"
        )
    }

    if (isTRUE(obsm_mapping) && !identical(spatialCoords_mapping, FALSE)) {
        spatial_key <- spatialCoords_mapping[[1]]
        rd_names <- SingleCellExperiment::reducedDimNames(spe)
        if (spatial_key %in% rd_names) {
            cli_abort(c(
                "Key {.val {spatial_key}} is reserved for spatial coordinates.",
                "i" = "Use a different {.arg obsm_mapping} target or set {.arg spatialCoords_mapping = FALSE}."
            ))
        }
    }

    if (isTRUE(uns_mapping) && isTRUE(imgData_mapping)) {
        meta_names <- names(S4Vectors::metadata(spe))
        reserved_uns <- c("spatial", "anndataR_spatialexperiment")
        if (any(reserved_uns %in% meta_names)) {
            cli_abort(c(
                "Reserved metadata keys detected in {.code metadata(spe)}.",
                "i" = "The following keys are reserved for SpatialExperiment export: {.val {reserved_uns}}"
            ))
        }
    }

    adata <- from_SingleCellExperiment(
        sce=spe,
        x_mapping=x_mapping,
        layers_mapping=layers_mapping,
        obs_mapping=obs_mapping,
        var_mapping=var_mapping,
        obsm_mapping=obsm_mapping,
        varm_mapping=varm_mapping,
        obsp_mapping=obsp_mapping,
        varp_mapping=varp_mapping,
        uns_mapping=uns_mapping,
        output_class=output_class,
        ...
    )

    if (!identical(spatialCoords_mapping, FALSE)) {
        .from_SPE_process_spatialCoords(
            adata=adata,
            spe=spe,
            spatial_key=spatialCoords_mapping[[1]]
        )
    }

    if (isTRUE(imgData_mapping)) {
        .from_SPE_process_imgData(
            adata=adata,
            spe=spe,
            library_id=library_id,
            image_mode=image_mode
        )
    }

    .from_SPE_store_roundtrip_metadata(
        adata=adata,
        spe=spe,
        spatialCoords_mapping=spatialCoords_mapping,
        library_id=library_id,
        image_mode=image_mode
    )

    adata
}


.from_SPE_process_spatialCoords <- function(adata, spe, spatial_key="spatial") {
    coords <- SpatialExperiment::spatialCoords(spe)

    if (is.null(coords) || ncol(coords) == 0L || nrow(coords) == 0L) {
        return(invisible())
    }

    coords <- as.matrix(coords)

    if (!is.numeric(coords)) {
        cli_abort(
            "{.code spatialCoords(spe)} must be numeric to be exported to {.code adata$obsm[[{spatial_key}]]}"
        )
    }

    rownames(coords) <- colnames(spe)
    adata$obsm[[spatial_key]] <- coords
    invisible()
}


.from_SPE_guess_library_ids <- function(spe, library_id="auto") {
    cd <- SummarizedExperiment::colData(spe)

    if (identical(library_id, "sample_id")) {
        if (!("sample_id" %in% names(cd))) {
            cli_abort("Column {.val sample_id} not found in {.code colData(spe)}")
        }
        return(as.character(cd[["sample_id"]]))
    }

    if (identical(library_id, "auto")) {
        if ("sample_id" %in% names(cd)) {
            return(as.character(cd[["sample_id"]]))
        }

        img <- tryCatch(SpatialExperiment::imgData(spe), error=function(e) NULL)
        if (!is.null(img) && nrow(img) > 0L && "sample_id" %in% names(img)) {
            vals <- unique(as.character(img[["sample_id"]]))
            if (length(vals) == 1L) {
                return(rep(vals, ncol(spe)))
            }
        }

        return(rep("library1", ncol(spe)))
    }

    if (length(library_id) == 1L && is.character(library_id)) {
        if (!(library_id %in% names(cd))) {
            cli_abort(
                "Column {.val {library_id}} not found in {.code colData(spe)}"
            )
        }
        return(as.character(cd[[library_id]]))
    }

    cli_abort(
        "{.arg library_id} must be one of {.val auto}, {.val sample_id}, or a column name in {.code colData(spe)}"
    )
}


.from_SPE_get_imgData <- function(spe) {
    tryCatch(
        SpatialExperiment::imgData(spe),
        error=function(e) NULL
    )
}


.from_SPE_image_to_uns <- function(x, image_mode) {
    if (identical(image_mode, "none")) {
        return(NULL)
    }

    if (identical(image_mode, "reference")) {
        if (is.character(x) && length(x) == 1L) {
            return(x)
        }
        return(NULL)
    }

    x
}


.from_SPE_split_img_row <- function(row_df) {
    nms <- names(row_df)

    image_col <- intersect(
        c("data", "image", "img", "raster", "loaded"),
        nms
    )
    image_col <- if (length(image_col) > 0L) image_col[[1]] else NULL

    image_id_col <- intersect(
        c("image_id", "imageID", "image", "type", "resolution"),
        nms
    )
    image_id_col <- if (length(image_id_col) > 0L) image_id_col[[1]] else NULL

    sample_col <- intersect(
        c("sample_id", "sample_idImg", "library_id"),
        nms
    )
    sample_col <- if (length(sample_col) > 0L) sample_col[[1]] else NULL

    scale_names <- c(
        "spot_diameter_fullres",
        "fiducial_diameter_fullres",
        "tissue_hires_scalef",
        "tissue_lowres_scalef"
    )

    list(
        image_col=image_col,
        image_id_col=image_id_col,
        sample_col=sample_col,
        scale_cols=intersect(scale_names, nms)
    )
}


.from_SPE_process_imgData <- function(
    adata,
    spe,
    library_id="auto",
    image_mode="embed"
) {
    img <- .from_SPE_get_imgData(spe)

    if (is.null(img) || nrow(img) == 0L) {
        return(invisible())
    }

    img <- as.data.frame(img, stringsAsFactors=FALSE)

    if (is.null(adata$uns)) {
        adata$uns <- list()
    }
    if (is.null(adata$uns[["spatial"]])) {
        adata$uns[["spatial"]] <- list()
    }

    row_info <- .from_SPE_split_img_row(img)

    for (i in seq_len(nrow(img))) {
        row_i <- img[i, , drop=FALSE]

        lib <- if (!is.null(row_info$sample_col)) {
            as.character(row_i[[row_info$sample_col]])
        } else {
            unique(.from_SPE_guess_library_ids(spe, library_id))[[1]]
        }

        if (is.null(adata$uns[["spatial"]][[lib]])) {
            adata$uns[["spatial"]][[lib]] <- list(
                images=list(),
                scalefactors=list(),
                metadata=list()
            )
        }

        image_id <- if (!is.null(row_info$image_id_col)) {
            as.character(row_i[[row_info$image_id_col]])
        } else {
            paste0("image_", i)
        }

        if (!is.null(row_info$image_col)) {
            image_payload <- .from_SPE_image_to_uns(
                row_i[[row_info$image_col]][[1]],
                image_mode=image_mode
            )

            if (!is.null(image_payload)) {
                adata$uns[["spatial"]][[lib]][["images"]][[image_id]] <- image_payload
            }
        }

        if (length(row_info$scale_cols) > 0L) {
            for (nm in row_info$scale_cols) {
                adata$uns[["spatial"]][[lib]][["scalefactors"]][[nm]] <- row_i[[nm]]
            }
        }

        meta_cols <- setdiff(
            names(row_i),
            c(
                row_info$image_col,
                row_info$image_id_col,
                row_info$sample_col,
                row_info$scale_cols
            )
        )

        if (length(meta_cols) > 0L) {
            adata$uns[["spatial"]][[lib]][["metadata"]][[image_id]] <- as.list(row_i[meta_cols])
        }
    }

    invisible()
}


.from_SPE_store_roundtrip_metadata <- function(
    adata,
    spe,
    spatialCoords_mapping="spatial",
    library_id="auto",
    image_mode="embed"
) {
    if (is.null(adata$uns)) {
        adata$uns <- list()
    }

    coords <- tryCatch(
        SpatialExperiment::spatialCoords(spe),
        error=function(e) NULL
    )

    coord_names <- if (!is.null(coords)) colnames(coords) else NULL

    adata$uns[["anndataR_spatialexperiment"]] <- list(
        schema_version="1.0.0",
        spatialCoords_colnames=coord_names,
        spatialCoords_key=if (identical(spatialCoords_mapping, FALSE)) NULL else spatialCoords_mapping,
        library_id_col=library_id,
        image_mode=image_mode
    )

    invisible()
}
