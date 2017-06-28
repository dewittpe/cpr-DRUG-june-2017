# influence-examples.R
#
# This file build graphics needed for the example influence weights.
library(cpr)
library(ggplot2)
library(dplyr)


bmat <- bsplines(x = seq(0, 6, length = 500), iknots = c(1.0, 1.5, 2.3, 4.0, 4.5))
theta <- c(1, 0, 3.5, 4.2, 3.7, -0.5, -0.7, 2, 1.5)
eg_cp <- cp(bmat, theta)

basis_plot <-
  plot(bmat, show_xi = TRUE, show_x = TRUE, color = TRUE, digits = 1) +
  theme(text = element_text(family = "Times", size = 10),
        legend.position = "none")

cp_plot0 <- basis_plot
cp_plot0[[which(names(cp_plot0) == "layers")]] <- NULL
cp_plot0 +
geom_line(mapping = aes(x = x, y = y), data = get_spline(eg_cp)$spline,
          inherit.aes = FALSE) +
theme(text = element_text(family = "Times", size = 10))

cp_plot0 <- basis_plot
cp_plot0[[which(names(cp_plot0) == "layers")]] <- NULL

eg_spline_and_cp_plot <-
  cp_plot0 +
    coord_cartesian(ylim = c(-1.3, 4.8)) +
  geom_line(mapping = aes(x = xi_star, y = theta), data = eg_cp$cp,
            inherit.aes = FALSE, linetype = 3) +
  geom_line(mapping = aes(x = x, y = y), data = get_spline(eg_cp)$spline,
            inherit.aes = FALSE) +
  geom_point(mapping = aes(x = xi_star, y = theta),
             data = eg_cp$cp, inherit.aes = FALSE) +
  theme(text = element_text(family = "Times", size = 10))
ggsave(plot = eg_spline_and_cp_plot, filename = "figure/original_spline.png", width = 4, height = 3, units = "in")



omit_xi <- influence_of(cp(bmat, theta), c(6, 8))

# influence_data: need to construct this data.frame for plotting.  The version
# of the influence plots used in section 5 of the manuscript is far simpler.
influence_data <-
  bind_rows(rename(eg_cp$cp, x = xi_star, y = theta),
            get_spline(eg_cp)$spline,
            rename(omit_xi$coarsened_cp[[1]]$cp, x = xi_star, y = theta) ,
            get_spline(omit_xi$coarsened_cp[[1]])$spline,
            rename(omit_xi$coarsened_cp[[2]]$cp, x = xi_star, y = theta),
            get_spline(omit_xi$coarsened_cp[[2]])$spline,
            rename(omit_xi$reinserted_cp[[1]]$cp, x = xi_star, y = theta),
            get_spline(omit_xi$reinserted_cp[[1]])$spline,
            rename(omit_xi$reinserted_cp[[2]]$cp, x = xi_star, y = theta),
            get_spline(omit_xi$reinserted_cp[[2]])$spline,
            .id = "object")

influence_data$CP <- NA
influence_data$CP <- gsub("^[1-2]$", "Original", influence_data$object)
influence_data$CP <- gsub("^[3-6]$", "Omit a knot", influence_data$CP)
influence_data$CP <- gsub("^([7-9]|10)$", "Reinsert the knot", influence_data$CP)
influence_data$CP <-
  factor(influence_data$CP,
         levels = c("Original", "Omit a knot", "Reinsert the knot"))

influence_data <-
  bind_rows(filter(influence_data, object %in% c(1:2, 3:4)),
            filter(influence_data, object %in% c(1:2, 3:4, 7:8)),
            filter(influence_data, object %in% c(1:2, 5:6)),
            filter(influence_data, object %in% c(1:2, 5:6, 9:10)),
            .id = "facets")

influence_data$facets <-
  factor(influence_data$facets,
         levels = 1:4,
         labels = c("Omitting ~ xi[6]", "Reinserting ~ xi[6]",
                    "Omitting ~ xi[8]", "Reinserting ~ xi[8]"))

influence_weight_plot <- plot(bmat, show_x = FALSE)
influence_weight_plot[[which(names(influence_weight_plot) == "layers")]] <- NULL

influence_weight_plot <-
  influence_weight_plot %+%
    influence_data +
    facet_wrap( ~ facets, labeller = label_parsed) +
    aes(x = x, y = y, color = CP, shape = CP, linetype = CP,
        group = object) +
    geom_point(data = filter(influence_data, object %in% seq(1, 9, by = 2))) +
    geom_line() +
    theme(axis.ticks.y     = element_blank(),
          axis.text.y      = element_blank(),
          panel.grid       = element_blank(),
          legend.position  = "bottom",
          legend.key.width = unit(0.5, "inch"),
          legend.title     = element_blank()) +
    coord_cartesian(ylim = c(-1.3, 4.8)) +
    theme(text = element_text(family = "Times", size = 10))



with(influence_data,  ftable(facets, CP, object))

eg_spline_and_cp_plot_omit_xi6 <-
  eg_spline_and_cp_plot +
    geom_line(data = filter(influence_data, object == "4"),
              mapping = aes(x = x, y = y),
              color = "red",
              inherit.aes = FALSE) +
    geom_point(data = filter(influence_data, object == "3"),
              mapping = aes(x = x, y = y),
              color = "red",
              inherit.aes = FALSE) +
    geom_line(data = filter(influence_data, object == "3"),
              mapping = aes(x = x, y = y),
              color = "red", linetype = 2,
              inherit.aes = FALSE)

eg_spline_and_cp_plot_omit_xi8 <-
  eg_spline_and_cp_plot +
    geom_line(data = filter(influence_data, object == "6"),
              mapping = aes(x = x, y = y),
              color = "red",
              inherit.aes = FALSE) +
    geom_point(data = filter(influence_data, object == "5"),
              mapping = aes(x = x, y = y),
              color = "red",
              inherit.aes = FALSE) +
    geom_line(data = filter(influence_data, object == "5"),
              mapping = aes(x = x, y = y),
              color = "red", linetype = 2,
              inherit.aes = FALSE)

eg_spline_and_cp_plot_reinsert_xi6 <-
  eg_spline_and_cp_plot_omit_xi6 +
    geom_line(data = filter(influence_data, object == "8"),
              mapping = aes(x = x, y = y),
              color = "blue",
              inherit.aes = FALSE) +
    geom_point(data = filter(influence_data, object == "7"),
              mapping = aes(x = x, y = y),
              color = "blue",
              inherit.aes = FALSE) +
    geom_line(data = filter(influence_data, object == "7"),
              mapping = aes(x = x, y = y),
              color = "blue", linetype = 4,
              inherit.aes = FALSE)

eg_spline_and_cp_plot_reinsert_xi8 <-
  eg_spline_and_cp_plot_omit_xi8 +
    geom_line(data = filter(influence_data, object == "10"),
              mapping = aes(x = x, y = y),
              color = "blue",
              inherit.aes = FALSE) +
    geom_point(data = filter(influence_data, object == "9"),
              mapping = aes(x = x, y = y),
              color = "blue",
              inherit.aes = FALSE) +
    geom_line(data = filter(influence_data, object == "9"),
              mapping = aes(x = x, y = y),
              color = "blue", linetype = 4,
              inherit.aes = FALSE)

ggsave(plot = eg_spline_and_cp_plot_omit_xi6,
       file = "figure/eg_spline_and_cp_plot_omit_xi6.png",
       width = 4, height = 3)
ggsave(plot = eg_spline_and_cp_plot_reinsert_xi6,
       file = "figure/eg_spline_and_cp_plot_reinsert_xi6.png",
       width = 4, height = 3)
ggsave(plot = eg_spline_and_cp_plot_omit_xi8,
       file = "figure/eg_spline_and_cp_plot_omit_xi8.png",
       width = 4, height = 3)
ggsave(plot = eg_spline_and_cp_plot_reinsert_xi8,
       file = "figure/eg_spline_and_cp_plot_reinsert_xi8.png",
       width = 4, height = 3)
