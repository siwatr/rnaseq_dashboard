# Shared grouped colour/group-by choices helper (R/utils_group_choices.R).

test_that("group_field_choices orders General -> This session -> Data metadata", {
  ch <- group_field_choices(c("condition", "batch"),
                            session_items = c("Suggested removal" = "__removal__"),
                            none = TRUE)
  expect_equal(names(ch), c("General", "This session", "Data metadata"))
  expect_equal(unname(ch$General), "__none__")
  expect_equal(unname(ch$`This session`), "__removal__")
  expect_equal(unname(ch$`Data metadata`), c("condition", "batch"))
})

test_that("group_field_choices omits General when none = FALSE", {
  ch <- group_field_choices(c("condition"),
                            session_items = c("In removal pool" = "__pool__"),
                            none = FALSE)
  expect_equal(names(ch), c("This session", "Data metadata"))
})

test_that("group_field_choices inserts a Spike-in optgroup between session + metadata", {
  ch <- group_field_choices(c("condition"),
                            session_items = c("Suggested removal" = "__removal__"),
                            none = TRUE,
                            spike_items = c("Dose-response slope" = "__spike__slope"))
  expect_equal(names(ch), c("General", "This session", "Spike-in", "Data metadata"))
  expect_equal(unname(ch$`Spike-in`), "__spike__slope")
})

test_that("group_field_choices flattens a lone Data-metadata group", {
  ch <- group_field_choices(c("condition", "batch"), session_items = NULL, none = FALSE)
  expect_false(is.list(ch))                          # bare named vector, no optgroup
  expect_equal(unname(ch), c("condition", "batch"))
  expect_equal(names(ch), c("condition", "batch"))
})

test_that("group_field_choices drops empty optgroups", {
  ch <- group_field_choices(character(0), session_items = NULL, none = TRUE)
  expect_equal(names(ch), "General")                 # no colData, no session
})
