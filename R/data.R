#' Trap positions from Bouyer et al. (Table 2)
#'
#' Coordinates (x, y) in metres of 21 mosquito traps used in the reference
#' study. The origin (0, 0) corresponds to the sterile mosquito release point.
#'
#' @format A numeric matrix with 21 rows and 2 columns: `x_coord` and `y_coord`.
#' @source Bouyer et al., Table 2.
"mosquito_trap_positions"


#' Daily mosquito captures from Bouyer et al. (Table 1)
#'
#' Number of mosquitoes captured per trap per day over 20 days, derived from
#' the cumulative counts in the reference article (daily = cumul[t] - cumul[t-1]).
#'
#' @format An integer matrix with 20 rows (days) and 21 columns (traps).
#' @source Bouyer et al., Table 1.
"mosquito_captures"
