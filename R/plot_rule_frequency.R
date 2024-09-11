#' @title
#' Plot frequency of rule discovery
#'
#' @description
#' Generates a plot showing the number of times each rule (before filtering)
#' was discovered.
#'
#' @param object A `cre` object.
#' @param max_rules Maximum number of decision rules to plot.
#'
#' @return
#' A plot showing the frequency of rule discovery.
#'

plot_rule_frequency <- function(object, max_rules = 10, ){

  `%>%` <- magrittr::`%>%`

  rule_freq_df <- object[["rule_freq_list"]] %>% as.data.frame()

  message("Before filtering, ", nrow(rule_freq_df), " rules were discovered.")

  if(max_rules == 10){
    message("Plotting 10 rules by default. Specify `max_rules` to plot
more or fewer rules.")
  }

  if(nrow(rule_freq_df > max_rules)){
    rule_freq_df <- rule_freq_df[1:max_rules,]
  }

  # clean up the way the rules are written if categorical
  rule_freq_df$rules <- gsub("%in%", "in", rule_freq_df$rules)

  # remove quotes around strings
  rule_freq_df$rules <- gsub("'", "", rule_freq_df$rules)

  # # remove c()
  # rule_freq_df$rules <- gsub("c(", "", rule_freq_df$rules)
  # rule_freq_df$rules <- gsub(")", "", rule_freq_df$rules)

  # make rules factors and order by frequency
  rule_freq_df$rules <- factor(rule_freq_df$rules,
                               levels = rev(rule_freq_df$rules))

  # plot
  rule_freq_df %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = freq, y = rules)) +
    ggplot2::labs(x = "% of iterations selected", y = "") +
    ggplot2::xlim(c(0,100)) +
    ggplot2::theme_bw()

}
