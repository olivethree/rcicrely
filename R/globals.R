## Package-level declarations.
##
## Both of these are required when a package uses `[.data.table`
## syntax but does not `import(data.table)` in its NAMESPACE (we use
## `data.table::` prefixes instead). See CLAUDE.md sec.9.1.
##
## - `utils::globalVariables()` silences R CMD check NOTEs about
##   unbound sentinels used inside data.table calls.
## - `.datatable.aware <- TRUE` flips the runtime flag that makes
##   `cedta` recognise this package as data.table-aware, so
##   `[.data.table` dispatches correctly instead of falling through to
##   `[.data.frame`.

utils::globalVariables(c(".N", ".SD", "response"))

.datatable.aware <- TRUE
