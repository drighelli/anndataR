test_that("SpatialExperiment -> AnnData exports spatial coordinates", {
    skip_if_not_installed("SpatialExperiment")

    counts <- matrix(seq_len(20), nrow=5, ncol=4)
    colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))

    coords <- cbind(
        x=c(10, 20, 30, 40),
        y=c(1, 2, 3, 4)
    )

    spe <- SpatialExperiment::SpatialExperiment(
        assays=list(counts=counts),
        colData=S4Vectors::DataFrame(
            sample_id=c("s1", "s1", "s1", "s1"),
            row.names=colnames(counts)
        ),
        rowData=S4Vectors::DataFrame(
            feature_id=rownames(counts),
            row.names=rownames(counts)
        ),
        spatialCoords=coords
    )

    adata <- as_AnnData(spe)

    expect_true("spatial" %in% adata$obsm_keys())
    expect_equal(unname(as.matrix(adata$obsm[["spatial"]])), unname(coords))
    expect_true("anndataR_spatialexperiment" %in% adata$uns_keys())
})

