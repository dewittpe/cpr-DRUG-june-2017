library(cpr)
library(cprtesting)
library(animation)
library(ggplot2)

fdat <- progesterone1d(1, 1000)
fdat_init_cp <- cp(progesterone ~ bsplines(day, df = 54), data = fdat)
fdat_cpr_run <- cpr(fdat_init_cp)

des = "CPR Example Run with known Target Function"

# ?saveHTML
saveHTML({
  for (i in rev(seq_along(fdat_cpr_run))) {

    this_cp <- get_spline(fdat_cpr_run[[i]])

    print(
          plot(fdat_cpr_run[[i]], show_spline = FALSE, show_cp = FALSE) +
          coord_cartesian(ylim = c(-0.5, 4.0)) +

          geom_line(data = this_cp$spline, mapping = aes(x = x, y = y, linetype = "Spline")) +

          geom_point(data = this_cp$cp, mapping = aes(x = xi_star, y = theta, color = "Spline & CP")) +
          geom_line(data = this_cp$cp, mapping = aes(xi_star, y = theta, color = "Spline & CP", linetype = "CP")) +

          geom_line(data = fdat,
                    mapping = aes(x = day, y = progesterone, color = "Target", linetype = "Target"),
                    inherit.aes = FALSE) +
          scale_color_manual(name = "", values = c("black", "red")) + 
          scale_linetype_manual(name = "", values = c("dotted", "solid", "dashed")) +
          theme(legend.position = "bottom",
                legend.key.size = unit(0.5, "inches"))
          )
    ani.pause()
  } 
},
img.name = "fdat_cpr_run_plot",
imgdir = "fdat_cpr_run_dir",
htmlfile = "fdat_cpr_run.html", 
autobrowse = FALSE,
title = des,
description = des,
interval = 0.10,
navigator = FALSE,
verbose = FALSE,
loop = TRUE,
autoplay = FALSE,
ani.width=600,
ani.height=400)

# What would a forward step CP look like?
cp_frwd <- cp(progesterone ~ bsplines(day, df = 11), data = fdat)

g <-
  plot(cp_frwd, fdat_cpr_run[[8]], show_cp = FALSE, show_spline = TRUE, color = TRUE) +
  geom_line(data = fdat,
            mapping = aes(x = day, y = progesterone, color = "Target", linetype = "Target"),
            inherit.aes = FALSE) +
  scale_color_manual(name = "",
                     values = c("blue", "black", "red"),
                     label = c("Knots on Quantiles",
                               "CPR Selected",
                               "Target")) +
  scale_linetype(name = "",
                 label = c("Knots on Quantiles",
                           "CPR Selected",
                           "Target")) +
  theme(legend.position = "bottom")
ggsave(filename = "figure/ideal_cpr_v_frwd.png", plot = g,
       width = 6, height = 4, units = "in")



