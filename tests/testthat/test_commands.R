context("Commands Log Parsing")

log <- basejump::import("surecell_commands.log")



# getBarcodeCutoffFromCommands =================================================
test_that("getBarcodeCutoffFromCommands", {
    expect_identical(
        object = getBarcodeCutoffFromCommands(log),
        expected = 1000L
    )
})



# getLevelFromCommands =========================================================
test_that("getLevelFromCommands", {
    expect_identical(
        object = getLevelFromCommands(log),
        "genes"
    )
})



# getUMITypeFromCommands =======================================================
test_that("getUMITypeFromCommands", {
    expect_identical(
        object = getUMITypeFromCommands(log),
        expected = "surecell"
    )
})
