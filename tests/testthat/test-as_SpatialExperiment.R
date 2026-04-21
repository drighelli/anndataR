
test_that("AnnData -> SpatialExperiment imports spatial coordinates", {
    skip_if_not_installed("SpatialExperiment")

    adata <- AnnData(
        X=matrix(seq_len(20), nrow=4, ncol=5),
        obs=data.frame(
            sample_id=c("s1", "s1", "s1", "s1"),
            row.names=paste0("cell", 1:4)
        ),
        var=data.frame(
            row.names=paste0("gene", 1:5)
        ),
        obsm=list(
            spatial=cbind(
                x=c(10, 20, 30, 40),
                y=c(1, 2, 3, 4)
            )
        ),
        uns=list(
            anndataR_spatialexperiment=list(
                schema_version="1.0.0",
                spatialCoords_colnames=c("x", "y"),
                spatialCoords_key="spatial"
            )
        )
    )

    spe <- as_SpatialExperiment(adata)

    expect_s4_class(spe, "SpatialExperiment")
    expect_equal(
        unname(as.matrix(SpatialExperiment::spatialCoords(spe))),
        matrix(c(10, 1, 20, 2, 30, 3, 40, 4), ncol=2, byrow=TRUE)
    )
    expect_equal(colnames(SpatialExperiment::spatialCoords(spe)), c("x", "y"))
})
