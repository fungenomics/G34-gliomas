


#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {
  
  theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
      axis.ticks = element_line(colour = border_colour),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black", size = rel(1.2)),
      strip.text.y = element_text(colour = "black", size = rel(1.2)),
      title = element_text(size = rel(0.9)),
      axis.text = element_text(colour = "black", size = rel(0.8)),
      axis.title = element_text(colour = "black", size = rel(0.9)),
      legend.title = element_text(colour = "black", size = rel(0.9)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}


#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotateX <- function(angle = 90) {
  
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
  
}




#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
noLegend <- function() {
  
  theme(legend.position = "none")
  
}
