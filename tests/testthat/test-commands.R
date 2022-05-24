log <- import(file.path("cache", "surecell-commands.log"))

test_that("getBarcodeCutoffFromCommands", {
    expect_identical(
        object = getBarcodeCutoffFromCommands(log),
        expected = 1000L
    )
})

test_that("getLevelFromCommands", {
    expect_identical(
        object = getLevelFromCommands(log),
        "genes"
    )
})

test_that("getUMITypeFromCommands", {
    expect_identical(
        object = getUMITypeFromCommands(log),
        expected = "surecell"
    )
})
