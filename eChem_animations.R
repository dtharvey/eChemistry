#' The functions here are used to animate the results of simulating cyclic voltammetry, linear sweep voltammetry, chronoamperometry, and chronocoulometery experiments. Each function operates on an object created using cvSim, lsvSim, caSim, or ccSim.

#' voltgramGIF: function that creates an animated gif of a cyclic voltammogram or linear sweep voltammogram

voltgramGIF = function(filename, gif.name = "vgram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use vgramGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options  
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames  
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save the gif
  saveGIF({
    for (i in seq(1, length(filename$potential), pot_increment)) {
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options  
  ani.options(old.ani)
}

#' ampgramGIF: function that creates an animated gif  of a chronoamperogram

ampgramGIF = function(filename, gif.name = "cagram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cagramGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options  
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames  
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save the gif
  saveGIF({
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$time[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "time (s)", ylab = "current (µA)")
      grid()
    }}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options  
  ani.options(old.ani)
}

#' coulgramGIF: function that creates an animated gif of a chronocoulogram

coulgramGIF = function(filename, gif.name = "ccgram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cagramGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options  
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames  
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save the gif
  saveGIF({
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$time[1:i], y = filename$charge[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$charge), max(filename$charge)),
           xlab = "time (s)", ylab = "charge (µC)")
      grid()
    }}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options  
  ani.options(old.ani)
}

#' diffusionGIF: function that creates an animated gif of the diffusion profiles for a cyclic voltammetry, linear sweep voltammetry, or chronoamperometry experiment

diffusionGIF = function(filename, gif.name = "diffusion"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use diffusionGif.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options   
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save the gif
  saveGIF({
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
    }}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options   
  ani.options(old.ani)
}

#' cvGIF: function that creates an animated gif a cyclic voltammetry experiment, including both the cyclic voltammogram and the corresponding diffuions profiles

cvGIF = function(filename, gif.name = "cv"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cvGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save the gif  
  saveGIF({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }
    par(old.par)}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options 
  ani.options(old.ani)
}

#' lsvGIF: function that creates an animated gif of a linear sweep voltammetry experiment, including both the linear sweep voltammogram and the corresponding diffuions profiles

lsvGIF = function(filename, gif.name = "lsv"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use lsvGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save the gif  
  saveGIF({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }
    par(old.par)}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options 
  ani.options(old.ani)
}

# caGIF: function that creates an animated gif of a chronoamperometry experiment, including both the chronoamperogram and the corresponding diffuions profiles

caGIF = function(filename, gif.name = "ca"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use caGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save the gif  
  saveGIF({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$time[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "time (s)", ylab = "current (µA)")
      grid()
    }
    par(old.par)}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options 
  ani.options(old.ani)
}

# ccGIF: function that creates an animated gif of a chronocoulometry experiment, including both the chronocoulogram and the corresponding diffuions profiles

ccGIF = function(filename, gif.name = "cc"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use caGIF.")
  }
  library(animation)
  
#' save original ani.options and adjust the interval and loop options
  old.ani = ani.options(interval = 0.2, loop = 1)
  
#' set the increment for creating frames
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save the gif  
  saveGIF({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$time[1:i], y = filename$charge[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$charge), max(filename$charge)),
           xlab = "time (s)", ylab = "charge (µC)")
      grid()
    }
    par(old.par)}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options 
  ani.options(old.ani)
}

#' voltgramHTML: function that creates an HTML file with a movie of a cyclic voltammogram or a linear sweep voltammogram

voltgramHTML = function(filename, html.name = "voltgram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use vgramHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames  
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save html files 
  saveHTML({
    for (i in seq(1, length(filename$potential), pot_increment)) {
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }}, 
    img.name = paste0(html.name,"_plot"), 
    imgdir = paste0(html.name,"_dir"),
    htmlfile = paste0(html.name,".html"), 
    navigator = FALSE
  )
  
#' reset the original ani.options  
  ani.options(old.ani)
  
}

#' ampgramHTML: function that creates an HTML file with a movie of a chronoamperogram

ampgramHTML = function(filename, html.name = "ampgram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cagramHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames  
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save html files 
  saveHTML({
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$time[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "time (s)", ylab = "current (µA)")
      grid()
    }}, 
    img.name = paste0(html.name,"_plot"), 
    imgdir = paste0(html.name,"_dir"),
    htmlfile = paste0(html.name,".html"), 
    navigator = FALSE
  )
  
#' reset the original ani.options  
  ani.options(old.ani)
  
}

#' coulgramHTML: function that creates an HTML file with a movie of a chronocoulogram

coulgramHTML = function(filename, html.name = "coulgram"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cagramHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames  
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save html files 
  saveHTML({
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$time[1:i], y = filename$charge[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$charge), max(filename$charge)),
           xlab = "time (s)", ylab = "charge (µC)")
      grid()
    }}, 
    img.name = paste0(html.name,"_plot"), 
    imgdir = paste0(html.name,"_dir"),
    htmlfile = paste0(html.name,".html"), 
    navigator = FALSE
  )
  
#' reset the original ani.options  
  ani.options(old.ani)
  
}

#' diffusionHTML: function that creates an HTML file with a movie of a diffusion profile for a cyclic voltammetry, linear sweep voltammetry, or chronoamperometry experiment

diffusionHTML = function(filename, html.name = "diffusion"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use diffusionHTML.")
  }
  library(animation)
  
#' adjust ani.options
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save html files
  saveHTML({
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i, ], 
           col = "blue", type = "l", lwd = 3,
           ylim = c(0, 1050 * filename$conc.bulk),
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      lines(x = filename$distance, y = filename$reddata[i, ], 
        lwd = 3, col = "red")
      grid()
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
    }}, 
    img.name = paste0(html.name,"_plot"), 
    imgdir = paste0(html.name,"_dir"), 
    htmlfile = paste0(html.name,".html"), 
    navigator = FALSE
  )
  
#' reset the original ani.options    
  ani.options(old.ani)
  
}

#' cvHTML: function that creates an HTML movie of a cyclic voltammetry experiment, including both the cyclic voltammogram and the corresponding diffuions profiles

cvHTML = function(filename, html.name = "cv"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use cvHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save html files
  saveHTML({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }
    par(old.par)
  }, 
  img.name = paste0(html.name,"_plot"), 
  imgdir = paste0(html.name,"_dir"), 
  htmlfile = paste0(html.name,".html"), 
  navigator = FALSE
  )
  
#' reset the original ani.options 
  ani.options(old.ani)
  
}

#' lsvHTML: function that creates an  HTML movie of a linear sweep voltammetry experiment, including both the linear sweep voltammogram and the corresponding diffuions profiles

lsvHTML = function(filename, html.name = "lsv"){
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use lsvHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames
  pot_increment = round(length(filename$potential)/40, digits = 0)
  
#' create and save html files
  saveHTML({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), pot_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$potential[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(max(filename$potential), min(filename$potential)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "potential (V)", ylab = "current (µA)")
      grid()
    }
    par(old.par)
  }, 
  img.name = paste0(html.name,"_plot"), 
  imgdir = paste0(html.name,"_dir"), 
  htmlfile = paste0(html.name,".html"), 
  navigator = FALSE
  )
  
#' reset the original ani.options 
  ani.options(old.ani)
  
}

#' caHTML: function that creates an HTML movide of a chronoamperometry experiment, including both the chronoamperogram and the corresponding diffuions profiles

caHTML = function(filename, html.name = "ca") {
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use caHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save html files
  saveHTML({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$time[1:i], y = filename$current[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$current), max(filename$current)),
           xlab = "time (s)", ylab = "current (µA)")
      grid()
    }
    par(old.par)
  }, 
  img.name = paste0(html.name,"_plot"), 
  imgdir = paste0(html.name,"_dir"), 
  htmlfile = paste0(html.name,".html"), 
  navigator = FALSE
  )
  
#' reset the original ani.options 
  ani.options(old.ani)
  
}

#' ccHTML: function that creates an HTML movie of a chronocoulometry experiment, including both the chronocoulogram and the corresponding diffuions profiles

ccHTML = function(filename, html.name = "cc") {
  
  if (!requireNamespace("animation", quietly = TRUE)) {
    stop("You need to install the aniimation package to use caHTML.")
  }
  library(animation)
  
#' adjust ani.options  
  old.ani = ani.options(interval = 0.2, verbose = FALSE)
  
#' set increment between frames
  time_increment = round(length(filename$time)/40, digits = 0)
  
#' create and save html files
  saveHTML({
    old.par = par(mfrow = c(2, 1))
    for (i in seq(1, length(filename$time), time_increment)) {
      plot(x = filename$distance, y = filename$oxdata[i,  ],
           type = "l", lwd = 3, col = "blue", 
           ylim = c(0, 1050 * filename$conc.bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism != "E") {
        lines(x = filename$distance, 
              y = filename$chemdata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("Ox", "Red", "Chem"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = 0.05)
      } else {
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = 0.05)
      }
      plot(x = filename$time[1:i], y = filename$charge[1:i], 
           col = "blue", type = "l", lwd = 3,
           xlim = c(min(filename$time), max(filename$time)),
           ylim = c(min(filename$charge), max(filename$charge)),
           xlab = "time (s)", ylab = "charge (µC)")
      grid()
    }
    par(old.par)
  }, 
  img.name = paste0(html.name,"_plot"), 
  imgdir = paste0(html.name,"_dir"), 
  htmlfile = paste0(html.name,".html"), 
  navigator = FALSE
  )
  
#' reset the original ani.options 
  ani.options(old.ani)
  
}
