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

test_that("getUmiTypeFromCommands", {
    expect_identical(
        object = getUmiTypeFromCommands(log),
        expected = "surecell"
    )
})
