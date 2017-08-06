#' The functions here are used to plot the results of simulating cyclic voltammetry, linear sweep voltammetry, chronoamperometry, and chronocoulometery experiments. Each function operates on an object created using cvSim, lsvSim, caSim, or ccSim.

#' plotPotential: function that plots the applied potential as a function of time; filename is an object created using cycvolt, linsweep, or chronoamp; passing a character string to main_title adds a title to the plot

plotPotential = function(filename, main_title = NULL){
  plot(x = filename$time, y = filename$potential, lwd = 2, col = "blue",
       type = "l", xlab = "time (sec)", ylab = "potential (V)",
       main = main_title)
  grid()
}

#' plotDiffusion: function that plots a set of diffusion profiles (concentration as a function of distance from the electrode's surface) for an object created using cycvolt, linsweep, or chronoamp at the specified time; for an EC mechanism, the concentration of the inert product is calculated using a mass balance between the initial concentration of analyte and the combined concentrations of Ox, Red, and the product; function includes a default title that gives the time and the potential for the diffusion profile

plotDiffusion = function(filename, t = 1) {
  
  if (t < min(filename$time) | t > max(filename$time)) {
    stop(paste0("Time is limited to a value between ", min(filename$time), " s and ", max(filename$time), " s."))
  }
  
  index = which.min(abs(filename$time - t))
  
  ymax = filename$conc.bulk * 1000
  
  if (filename$mechanism == "E") {
    plot(x = filename$distance, y = filename$oxdata[index, ], 
         type = "l", lwd = 3, 
         col = "blue", ylim = c(0, ymax),
         xlab = "distance from electrode (cm)", 
         ylab = "concentration (mM)", 
         main = paste0(round(filename$time[index], digits = 3),
                       " sec & ", 
                       round(filename$potential[index], 
                             digits = 3), " V"))
    lines(x = filename$distance, y = filename$reddata[index, ], 
          lwd = 3, col = "red")
    legend(x = "right", legend = c("Ox", "Red"), 
           fill = c("blue", "red"), border = "white",
           bty = "n", inset = 0.05)
    
    grid()
  } else {
    plot(x = filename$distance, y = filename$oxdata[index, ], 
         type = "l", lwd = 3, 
         col = "blue", ylim = c(0, ymax),
         xlab = "distance from electrode (cm)", 
         ylab = "concentration (mM)", 
         main = paste0(round(filename$time[index], digits = 3),
                       " sec & ", 
                       round(filename$potential[index], 
                             digits = 3), " V"))
    lines(x = filename$distance, y = filename$reddata[index, ], 
          lwd = 3, col = "red")
    lines(x = filename$distance, 
          y = filename$chemdata[index, ], 
          lwd = 3, col = "green")
    legend(x = "right", legend = c("Ox", "Red", "Chem"), 
           fill = c("blue", "red", "green"), 
           bty = "n", inset = 0.05)
    grid()
  } 
}

#' plotGrid: function that plots eight diffusion profiles---at times that are 10%, 20%, 30%, 40%, 60%, 70%, 80%, and 90% of the total time for the experiment---around a central plot that shows the voltammogram or chronoamperogram from an object created useig cycvolt, linsweep, or chronoamp; all plots have default titles

plotGrid = function(filename) {
  
#' adjust the graphical parameters to create a 3 by 3 grid of graphical windows and to adjust the margins for each graphical window
  
  par(mfrow = c(3,3), mar = c(4.1, 4.1, 2.1, 2.1))
  
#' set the times for the diffusion profiles 
  
  t = round(c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9) * (filename$tunits) + 1, digits = 0)
  
  
# screen(1): diffusion profile at 20% of scan time
  
  plotDiffusion(filename, t = filename$time[t[2]])
  
# screen(2): diffusion profile at 30% of scan time
  
  plotDiffusion(filename, t = filename$time[t[3]])
  
# screen(3): diffusion profile at 40% of scan time
  
  plotDiffusion(filename, t = filename$time[t[4]])
  
# screen(4): diffusion profile at 10% of scan time
  
  plotDiffusion(filename, t = filename$time[t[1]])
  
# screen(5): voltammogram or chronoamperogram
  
  if (filename$expt == "CV" | filename$expt == "LSV") {
    plot(x = filename$potential, y = filename$current, 
         lwd = 2, col = "blue", type = "l", 
         xlab = "potential (V)", ylab = "current (µA)", 
         xlim = c(max(filename$potential), min(filename$potential)),
         ylim = c(1.1 * min(filename$current), 1.1 * max(filename$current)),
         main = "points show times")
    points(x = filename$potential[t + 1], 
           y = filename$current[t + 1], pch = 19, 
           col = "blue", cex = 1.5)
    for (i in 1:8) {
      if(t[i] < length(filename$time)/2){
        text(x = filename$potential[t[i] + 1], 
             y = filename$current[t[i] + 1], 
             labels = as.character(filename$time[t[i]]), 
             pos = 3, cex = 0.75)
      } else {
        text(x = filename$potential[t[i] + 1], 
             y = filename$current[t[i] + 1], 
             labels = as.character(filename$time[t[i]]), 
             pos = 1, cex = 0.75)
      }
    }
    grid()
  } else if (filename$expt == "CA") {
    plot(x = filename$time, y = filename$current, 
         lwd = 2, col = "blue", type = "l", 
         xlab = "timel (s)", ylab = "current (µA)", 
         main = "points show times")
    points(x = filename$time[t + 1], 
           y = filename$current[t + 1], 
           pch = 19, col = "blue", cex = 1.5)
    for (i in 1:8) {
      text(x = filename$time[t[i] + 1], 
           y = filename$current[t[i] + 1], 
           labels = as.character(filename$time[t[i]]), 
           pos = 3, cex = 0.75)
    }
    grid()
  } else {
    plot(x = filename$time, y = filename$charge, 
         lwd = 2, col = "blue", type = "l", 
         xlab = "timel (s)", ylab = "charge (µC)", 
         main = "points show times")
    points(x = filename$time[t + 1], 
           y = filename$charge[t + 1], 
           pch = 19, col = "blue", cex = 1.5)
    for (i in 1:8) {
      text(x = filename$time[t[i] + 1], 
           y = filename$charge[t[i] + 1], 
           labels = as.character(filename$time[t[i]]), 
           pos = 3, cex = 0.75)
    }
    grid()
  }
  
# screen(6): diffusion profiles at 60% of scan time
  
  plotDiffusion(filename, t = filename$time[t[5]])
  
# screen(7): diffusion profiles at 90% of scan time
  
  plotDiffusion(filename, t = filename$time[t[8]])
  
# screen(8): diffusion profiles at 80% of scan time
  
  plotDiffusion(filename, t = filename$time[t[7]])
  
# screen(9): diffusion profiles at 70% of scan time
  
  plotDiffusion(filename, t = filename$time[t[6]])
  
#' reestablish original graphical parameters
  
  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  
}

#' plotVoltgram: function that plots 1--5 voltammograms (cyclic or linear sweep) on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; line widths, line types, and line colors have default values that can be adjusted

plotVoltgram = function(filenames = list(file1, file2, file3, file, file5), 
                        legend_text = NULL, 
                        main_title = NULL,
                        line_widths = c(2, 2, 2, 2, 2), 
                        line_types = c(1, 2, 3, 4, 5),
                        line_colors = c("blue", "blue", "blue", 
                                     "blue", "blue")) {
  
  if (!requireNamespace("shape", quietly = TRUE)) {
    stop("You need to install the shape package to use plotVgram.")
  }
  
  library(shape) 
  
#' determine the number of voltammograms and verify that the number does not exceed the defined limit of five  
  numfiles = length(filenames)
  if (numfiles > 5){
    stop("Function is limited to plotting five voltammograms")
  }
  
#' determine the minimum and maximum units for the x-axis and the y-axis by examining the data for all voltammograms  
  xmin = rep(0, numfiles)
  xmax = rep(0, numfiles)
  ymin = rep(0, numfiles)
  ymax = rep(0, numfiles)
  
  for (i in 1:numfiles) {
    xmin[i] = min(filenames[[i]]$potential)
    xmax[i] = max(filenames[[i]]$potential)
    ymin[i] = min(filenames[[i]]$current)
    ymax[i] = max(filenames[[i]]$current)
  }
  
  xmin_id = which.min(xmin)
  xmax_id = which.max(xmax)
  ymin_id = which.min(ymin)
  ymax_id = which.max(ymax)
  
#' create the plot using lines for full data and points for reduced data
  if (filenames[[1]]$file_type == "full") {
    plot(x = filenames[[1]]$potential, 
         y = filenames[[1]]$current,
         xlim = c(xmax[xmax_id], xmin[xmin_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "l", lwd = line_widths[1] , lty = line_types[1], 
         col = line_colors[1],
         main = main_title,
         xlab = "potential (V)", ylab = "current (µA)"
    )
    Arrowhead(x0 = filenames[[1]]$potential[20], 
              y0 = filenames[[1]]$current[1],
              arr.type = "triangle", lcol = line_colors[1],
              angle = if(filenames[[1]]$direction == -1) {
                0
              } else {
                180
              })
  } else {
    plot(x = filenames[[1]]$potential, 
         y = filenames[[1]]$current,
         xlim = c(xmax[xmax_id], xmin[xmin_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "p", lwd = line_widths[1] , lty = line_types[1], 
         pch = 19, col = line_colors[1],
         xlab = "potential (V)", ylab = "current (µA)"
    )
  }
  
  grid()
  
  if (numfiles > 1) {
    for (i in 2:numfiles) {
      if (filenames[[i]]$file_type == "full") {
        lines(x = filenames[[i]]$potential, 
              y = filenames[[i]]$current, 
              lty = line_types[i], lwd = line_widths[i], 
              col = line_colors[i])
        Arrowhead(x0 = filenames[[i]]$potential[20], 
                  y0 = filenames[[i]]$current[1],
                  arr.type = "triangle", lcol = line_colors[i],
                  angle = if(filenames[[i]]$direction == -1) {
                    0
                  } else {
                    180
                  })
      } else {
        points(x = filenames[[i]]$potential, 
               y = filenames[[i]]$current, 
               pch = 19, col = line_colors[i])
      }
    }
  }
  if (is.null(legend_text) == FALSE) {
    legend(x = "topleft", legend = legend_text, lwd = line_widths, 
           lty = line_types, col = line_colors,bty = "n")
  }
}

#' plotAmpgram: function that plots 1--5 chronoamperograms on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; setting scale to a value less than 1 adjusts the y-axis limits so that the limits are not set by the current spikes; line widths, line types, and line colors have default values that can be adjusted

plotAmpgram = function(filenames = list(file1, file2, file3, file, file5),
                       scale = 1,
                       legend_text = NULL, 
                       main_title = NULL,
                       line_widths = c(2, 2, 2, 2, 2), 
                       line_types = c(1, 2, 3, 4, 5),
                       line_colors = c("blue", "blue", "blue", 
                                      "blue", "blue")) {
  
#' determine the number of chronoamperograms and verify that the number does not exceed the defined limit of five  
  
  numfiles = length(filenames)
  if (numfiles > 5){
    stop("Function is limited to plotting five chronoamperograms")
  }
  
#' determine the minimum and maximum units for the x-axis and the y-axis by examining the data for all chronoamperograms  
  
  xmin = rep(0, numfiles)
  xmax = rep(0, numfiles)
  ymin = rep(0, numfiles)
  ymax = rep(0, numfiles)
  
  for (i in 1:numfiles) {
    xmin[i] = min(filenames[[i]]$time)
    xmax[i] = max(filenames[[i]]$time)
    ymin[i] = min(filenames[[i]]$current) * scale
    ymax[i] = max(filenames[[i]]$current) * scale
  }
  
  xmin_id = which.min(xmin)
  xmax_id = which.max(xmax)
  ymin_id = which.min(ymin)
  ymax_id = which.max(ymax)
  
#' create the plot using lines for full data and points for reduced data
  
  if (filenames[[1]]$file_type == "full") {
    plot(x = filenames[[1]]$time, 
         y = filenames[[1]]$current,
         xlim = c( xmin[xmin_id], xmax[xmax_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "l", lwd = line_widths[1] , lty = line_types[1], 
         col = line_colors[1],
         main = main_title,
         xlab = "time (s)", ylab = "current (µA)"
    )
  } else {
    plot(x = filenames[[1]]$time, 
         y = filenames[[1]]$current,
         xlim = c( xmin[xmin_id], xmax[xmax_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "p", lwd = line_widths[1], lty = line_types[1], 
         pch = 19, col = line_colors[1],
         xlab = "time (s)", ylab = "current (µA)"
    )
  }
  
  grid()
  
  if (numfiles > 1) {
    for (i in 2:numfiles) {
      if (filenames[[i]]$file_type == "full") {
        lines(x = filenames[[i]]$time, 
              y = filenames[[i]]$current, 
              lty = line_types[i], lwd = line_widths[i], 
              col = line_colors[i])
      } else {
        points(x = filenames[[i]]$time, 
               y = filenames[[i]]$current, 
               pch = 19, col = line_colors[i])
      }
    }
  }
  if (is.null(legend_text) == FALSE) {
    legend(x = "topleft", legend = legend_text, lwd = line_widths, 
           lty = line_types, col = line_colors,bty = "n")
  }
}

#' plotCoulgram: function that plots 1--5 chronocoulomograms on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; setting scale to a value less than 1 adjusts the y-axis limits so that the limits are not set by the current spikes; line widths, line types, and line colors have default values that can be adjusted

plotCoulgram = function(filenames = list(file1, file2, file3, file, file5),
                        scale = 1,
                        legend_text = NULL, 
                        main_title = NULL,
                        line_widths = c(2, 2, 2, 2, 2), 
                        line_types = c(1, 2, 3, 4, 5),
                        line_colors = c("blue", "blue", "blue", 
                                      "blue", "blue")) {
  
#' determine the number of chronocoulograms and verify that the number does not exceed the defined limit of five  
  
  numfiles = length(filenames)
  if (numfiles > 5){
    stop("Function is limited to plotting five chronocouloograms")
  }
  
#' determine the minimum and maximum units for the x-axis and the y-axis by examining the data for all chronoamperograms  
  
  xmin = rep(0, numfiles)
  xmax = rep(0, numfiles)
  ymin = rep(0, numfiles)
  ymax = rep(0, numfiles)
  
  for (i in 1:numfiles) {
    xmin[i] = min(filenames[[i]]$time)
    xmax[i] = max(filenames[[i]]$time)
    ymin[i] = min(filenames[[i]]$charge) * scale
    ymax[i] = max(filenames[[i]]$charge) * scale
  }
  
  xmin_id = which.min(xmin)
  xmax_id = which.max(xmax)
  ymin_id = which.min(ymin)
  ymax_id = which.max(ymax)
  
#' create the plot using lines for full data and points for reduced data
  
  if (filenames[[1]]$file_type == "full") {
    plot(x = filenames[[1]]$time, 
         y = filenames[[1]]$charge,
         xlim = c( xmin[xmin_id], xmax[xmax_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "l", lwd = line_widths[1] , lty = line_types[1], 
         col = line_colors[1],
         main = main_title,
         xlab = "time (s)", ylab = "charge (µC)"
    )
  } else {
    plot(x = filenames[[1]]$time, 
         y = filenames[[1]]$charge,
         xlim = c( xmin[xmin_id], xmax[xmax_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "p", lwd = line_widths[1] , lty = line_types[1], 
         pch = 19, col = line_colors[1],
         xlab = "time (s)", ylab = "charge (µC)"
    )
  }
  
  grid()
  
  if (numfiles > 1) {
    for (i in 2:numfiles) {
      if (filenames[[i]]$file_type == "full") {
        lines(x = filenames[[i]]$time, 
              y = filenames[[i]]$charge, 
              lty = line_types[i], lwd = line_widths[i], 
              col = line_colors[i])
      } else {
        points(x = filenames[[i]]$time, 
               y = filenames[[i]]$charge, 
               pch = 19, col = line_colors[i])
      }
    }
  }
  if (is.null(legend_text) == FALSE) {
    legend(x = "topleft", legend = legend_text, lwd = line_widths, 
           lty = line_types, col = line_colors,bty = "n")
  }
}

#' annotateCV: function that plots a cyclic voltammogram and annotates it with values for Epc, Epa, delta E, ip,c, ip,a, and the peak current ratio; the percentage of points used to determine the peak currents is set using cathodic.per and anodic.per; to be measurable, the return peak current must satisfy a minimum threshold value and must be greater than the predicted baseline current 

annotateCV = function(filename, 
                      forward.per = 5, reverse.per = 5, 
                      threshold = 0.05) {
  
  if (filename$expt != "CV") {
    stop("This file is not from a cyclic voltammetry simulation.")
  }
  
  if (!requireNamespace("shape", quietly = TRUE)) {
    stop("You need to install the shape package to use annotateCV.")
  }
  library(shape)
  
#' plot the cyclic voltammogram
  
  plot(x = filename$potential, 
       y = filename$current, type = "l", 
       xlab = "potential (V)", ylab = "current (µA)", 
       xlim = c(max(filename$potential), min(filename$potential)))
  Arrowhead(x0 = filename$potential[20], 
            y0 = filename$current[1],
            arr.type = "triangle", 
            angle = if(filename$direction == -1) {
              0
            } else {
              180
            })
  grid()
  
#' identify Epc and Epa
  
  epc.id = which.max(filename$current)
  epa.id = which.min(filename$current)
  epc = filename$potential[epc.id]
  epa = filename$potential[epa.id]
  
#' create linear models to predict baseline currents using parameters forward.per and reverse.per; lm1 models the baseline for the foward reaction and uses the first forward.per% of the data points to model the baseline and a first-order relationship between current and potential; lm2 models the baseline for the reverse reaction and uses the first reverse.per% of the data points beginning with the switching potential to predict the current decay in the absence of a change in potential using a model in which the current decays as function of t^(-0.5)
  
  n.points1 = (forward.per/100) * filename$tunits
  n.points2 = (reverse.per/100) * filename$tunits
  lm1 = lm(filename$current[1:n.points1] ~ filename$potential[1:n.points1])
  lm2 = lm(filename$current[(filename$tunits/2+1):(filename$tunits/2+n.points2 + 1)] ~ I(filename$time[(filename$tunits/2 + 1):(filename$tunits/2+n.points2 + 1)]^(-0.5)))
  
#' add a dashed line to show the baseline for the forward reaction and add arrow to show the peak current and peak potential
  
  abline(lm1, lwd = 1, lty = 2, col = "black")
  if (filename$direction == -1) {
    arrows(x0 = epc, y0 = lm1$coefficients[1] + lm1$coefficients[2] * filename$potential[epc.id], 
           x1 = epc, y1 = filename$current[epc.id], 
           code = 3, length = 0.1, angle = 15)
  } else {
    arrows(x0 = epa, y0 = lm1$coefficients[1] + lm1$coefficients[2] * filename$potential[epc.id], 
           x1 = epa, y1 = filename$current[epa.id], 
           code = 3, length = 0.1, angle = 15)
  }
  
#' add a dashed line to show the baseline for the reverse reaction and, if threshold conditions are met, add arrow to show the peak current and peak potential
  
  y = lm2$coefficients[1] + lm2$coefficients[2] * filename$time[(filename$tunits/2 + 1):(filename$tunits+1)]^(-0.5)
  lines(filename$potential[(filename$tunits/2+1):(filename$tunits+1)], y,
        lwd = 1, lty = 2, col = "black")
  
#' set flags to control the plotting of annotations; flags for the cathodic and the anodic potentials determine whether potentials are included; flags for the cathodic and the anodic current determine whether currents are included; both sets of flags determine whether arrows are drawn
  
  flag.cpot = FALSE
  flag.ccur = FALSE
  flag.apot = FALSE
  flag.acur = FALSE
  
  if (abs(filename$current[epc.id]) < threshold) {
    flag.cpot = TRUE
    flag.ccur = TRUE
  }
  if (abs(filename$current[epa.id]) < threshold) {
    flag.apot = TRUE
    flag.acur = TRUE
  }
  
  if (filename$direction == -1 & flag.apot == FALSE) {
    if (filename$current[epa.id] - (y[epa.id - filename$tunits/2 + 1]) > 0) {
      flag.acur = TRUE
    }
  }
  if (filename$direction == 1 & flag.cpot == FALSE) {
    if (filename$current[epc.id] - (y[epc.id - filename$tunits/2 + 1]) < 0) {
      flag.ccur = TRUE
    }
  }
  
  if (filename$direction == -1) {
    if (flag.acur == FALSE & flag.apot == FALSE) {
      arrows(x0 = epa, y0 = y[epa.id - filename$tunits/2 + 1],
             x1 = epa, y1 = filename$current[epa.id],
             code = 3, length = 0.1, angle = 15)
    }
  } else {
    if (flag.ccur == FALSE & flag.cpot == FALSE) {
      arrows(x0 = epc, y0 = y[epc.id - filename$tunits/2 + 1],
             x1 = epc, y1 = filename$current[epc.id],
             code = 3, length = 0.1, angle = 15)
    }
  }
  
#' calculate ipc and ipa
  
  if (filename$direction == -1) {
    ipc = filename$current[epc.id] - (lm1$coefficients[1] + lm1$coefficients[2] * filename$potential[epc.id])
    ipa = filename$current[epa.id] - y[epa.id - filename$tunits/2+1]
  } else {
    ipa = (filename$current[epa.id] - (lm1$coefficients[1] + lm1$coefficients[2] * filename$potential[epa.id]))
    ipc = (filename$current[epc.id] - y[epc.id - filename$tunits/2+1])
  }
  
#' add annotations to plot for peak potentials, peak currents, delta E and peak current ratio; values are not shown for the reverse reaction if threshold values are not met
  
  delta.y = max(filename$current) - min(filename$current)
  
  #' annotations for cathodic peak potential and current
  
  if (flag.cpot == FALSE) {
    text(x = max(filename$potential), 
         y = max(filename$current),
         substitute(paste(E[pc], ": ", epc, " V"), list(epc = noquote(formatC(epc, format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    text(x = max(filename$potential), 
         y = max(filename$current),
         substitute(paste(E[pc], ": not measurable")), 
         adj = c(0, NA), cex = 0.80)
  }
  
  if (flag.ccur == FALSE) {
    text(x = min(filename$potential),
         y = max(filename$current) - 0.80 * delta.y,
         substitute(paste(i[pc], ": ", ipc, " µA"), list(ipc = noquote(formatC(ipc, format = "f", digits = 2)))), 
         adj = c(1, NA), cex = 0.80)
  } else {
    text(x = min(filename$potential),
         y = max(filename$current) - 0.80 * delta.y,
         substitute(paste(i[pc], ": not measurable")), 
         adj = c(1, NA), cex = 0.80)
  }
  
#' add annotations for anodic peak potential and current
  
  if (flag.apot == FALSE) {
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.05 * delta.y,
         substitute(paste(E[pa], ": ", epa, " V"), list(epa = noquote(formatC(epa, format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.05 * delta.y,
         substitute(paste(E[pa], ": not measurable")), 
         adj = c(0, NA), cex = 0.80)
  }
  
  if(flag.acur == FALSE) {
    text(x = min(filename$potential),
         y = max(filename$current) - 0.85 * delta.y,
         substitute(paste(i[pa], ": ", ipa, " µA"), list(ipa = noquote(formatC(ipa, format = "f", digits = 2)))), 
         adj = c(1, NA), cex = 0.80)
  } else {
    text(x = min(filename$potential),
         y = max(filename$current) - 0.85 * delta.y,
         substitute(paste(i[pa], ": not measurable")), 
         adj = c(1, NA), cex = 0.80)
  }
  
#' add annotations for delta E, Eavg, and current ratio
  
  if (flag.cpot == FALSE & flag.apot == FALSE) {
    text(x = max(filename$potential),
         y = max(filename$current) - 0.1 * delta.y,
         substitute(paste(Delta, "E: ", deltae, " V"),
                    list(deltae = noquote(formatC(epa - epc, format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = max(filename$potential),
         y = max(filename$current) - 0.15 * delta.y,
         substitute(paste(E[avg], ": ", eavg, " V"),
                    list(eavg = noquote(formatC(0.5 * (epa + epc), format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    text(x = max(filename$potential),
         y = max(filename$current) - 0.1 * delta.y,
         substitute(paste(Delta, "E: not measurable")),
         adj = c(0, NA), cex = 0.80)
    text(x = max(filename$potential),
         y = max(filename$current) - 0.15 * delta.y,
         substitute(paste(E[avg], ": not measurable")),
         adj = c(0, NA), cex = 0.80)
  }
  
  if (flag.ccur == FALSE & flag.acur == FALSE) {
    if (filename$direction == 1) {
      text(x = min(filename$potential),
           y = max(filename$current) - 0.90 * delta.y,
           substitute(paste("|", i[pc]/i[pa], "|: ", ratio), list(ratio = noquote(formatC(abs(ipc/ipa), format = "f", digits = 2)))), 
           adj = c(1, NA), cex = 0.80)
    } else {
      text(x = min(filename$potential),
           y = max(filename$current) - 0.90 * delta.y,
           substitute(paste("|", i[pa]/i[pc], "|: ", ratio), list(ratio = noquote(formatC(abs(ipa/ipc), format = "f", digits = 2)))), 
           adj = c(1, NA), cex = 0.80)
    }
  } else {
    text(x = min(filename$potential),
         y = max(filename$current) - 0.90 * delta.y,
         substitute(paste("|", i[pc]/i[pa], "|: not measurable")), 
         adj = c(1, NA), cex = 0.80)
  }
  
}

#' annotateLSV: function that plots a linear-sweep voltammogram and annotates it with values for either Epc and ip,c, or for Epa and ip,a,; the percentage of points used to determine the peak current is set using potential.per; see annotateCV for additional comments on the code

annotateLSV = function(filename, potential.per = 5) {
  
  if (filename$expt != "LSV") {
    stop("This file is not from a linear sweep voltammetry simulation.")
  }
  
  plot(x = filename$potential, y = filename$current, 
       type = "l", 
       xlab = "potential (V)", ylab = "current (µA)", 
       xlim = c(max(filename$potential), min(filename$potential)))
  grid()
  
  if (filename$direction == -1) {
    epc.id = which.max(filename$current)
    epc = filename$potential[epc.id]
  } else {
    epa.id = which.min(filename$current)
    epa = filename$potential[epa.id]
  }
  
  n.points = (potential.per/100) * filename$tunits
  lm = lm(filename$current[1:n.points] ~ filename$potential[1:n.points])
  abline(lm, lwd = 1, lty = 2, col = "black") 
  
  if (filename$direction == -1) {
    arrows(x0 = epc, y0 = lm$coefficients[1] + lm$coefficients[2] * filename$potential[epc.id], x1 = epc, y1 = filename$current[epc.id], code = 3, length = 0.1, angle = 15)
  } else {
    arrows(x0 = epa, y0 = lm$coefficients[1] + lm$coefficients[2] * filename$potential[epa.id], x1 = epa, y1 = filename$current[epa.id], code = 3, length = 0.1, angle = 15)
  }
  
  if (filename$direction == -1) {
    ipc = filename$current[epc.id] - (lm$coefficients[1] + lm$coefficients[2] * filename$potential[epc.id])
  } else {
    ipa = -(filename$current[epa.id] - (lm$coefficients[1] + lm$coefficients[2] * filename$potential[epa.id]))
  }
  
  delta.y = max(filename$current) - min(filename$current)
  if (filename$direction == -1) {
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.05 * delta.y,
         substitute(paste(E[pc], ": ", epc, " V"), list(epc = noquote(formatC(epc, format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.1 * delta.y,
         substitute(paste(i[pc], ": ", ipc, " µA"), list(ipc = noquote(formatC(ipc, format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.05 * delta.y,
         substitute(paste(E[pa], ": ", epa, " V"), list(epa = noquote(formatC(epa, format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = max(filename$potential), 
         y = max(filename$current) - 0.1 * delta.y,
         substitute(paste(i[pa], ": ", ipa, " µA"), list(ipa = noquote(formatC(ipa, format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
  }
}

#' annotateCA: function that plots a chronoamperogam and annotates it with either the current for a single pulse experiment, or with the current following the forward and the reverse pulse, and the current ratio for a double pulse experiment; time.delay is the time after a pulse at which the current is measured, which default to the length of the pulse if a value is not provided; increasing the scale.factor narrows the range of values on the y-axis

annotateCA = function(filename, time.delay, scale.factor = 1) {
  
  if (missing(time.delay)) {
    time.delay = filename$time_pulse2 - filename$time_pulse1
  }
  
  if (filename$expt != "CA") {
    stop("This file is not from a chronoamperometry simulation.")
  }
  
  plot(filename$time, filename$current, type = "l",
       xlab = "time (s)", ylab = "current (µA)",
       ylim = c(min(filename$current)/scale.factor, max(filename$current)/scale.factor))
  grid()
  
  delta.y = max(filename$current)/scale.factor - min(filename$current)/scale.factor
  
  abline(h = 0, lty = 2)
  if (filename$pulses == 1) {
    index = filename$tunits * (filename$time_pulse1 + time.delay)/filename$time_end
    arrows(x0 = filename$time[index], y0 = 0,
           x1 = filename$time[index], y1 = filename$current[index],
           code = 3, length = 0.1, angle = 15)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$current)/scale.factor - 0.10 * delta.y,
         substitute(paste("i: ", i, " µA"), list(i = noquote(formatC(filename$current[index], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    index1 = filename$tunits * (filename$time_pulse1 + time.delay)/filename$time_end
    index2 = filename$tunits * (filename$time_pulse2 + time.delay)/filename$time_end
    arrows(x0 = filename$time[index1], y0 = 0,
           x1 = filename$time[index1], y1 = filename$current[index1],
           code = 3, length = 0.1, angle = 15)
    arrows(x0 = filename$time[index2], y0 = 0,
           x1 = filename$time[index2], y1 = filename$current[index2],
           code = 3, length = 0.1, angle = 15)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$current)/scale.factor - 0.10 * delta.y,
         substitute(paste(i[f], ": ", ifor, " µA"), list(ifor = noquote(formatC(filename$current[index1], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$current)/scale.factor - 0.175 * delta.y,
         substitute(paste(i[r], ": ", irev, " µA"), list(irev = noquote(formatC(filename$current[index2], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$current)/scale.factor - 0.25 * delta.y,
         substitute(paste("|",i[r]/i[f], "|: ", ratio), list(ratio = noquote(formatC(abs(filename$current[index2]/filename$current[index1]), format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
  }
}

#' annotateCC: function that plots a chronocoulogam and annotates it with either the charge for a single pulse experiment, or with the charge following the forward and the reverse pulse, and the charge ratio for a double pulse experiment; time.delay is the time after a pulse at which the current is measured, which defaults to the length of the pulse if a value is not provided; increasing the scale.factor narrows the range of values on the y-axis

annotateCC = function(filename, time.delay, scale.factor = 1) {
  
  if (missing(time.delay)) {
    time.delay = filename$time_pulse2 - filename$time_pulse1
  }
  
  if (filename$expt != "CC") {
    stop("This file is not from a chronocoulometry simulation.")
  }
  
  plot(filename$time, filename$charge, type = "l",
       xlab = "time (s)", ylab = "charge (µC)",
       ylim = c(min(filename$charge)/scale.factor, max(filename$charge)/scale.factor))
  grid()
  
  delta.y = max(filename$charge)/scale.factor - min(filename$charge)/scale.factor
  
  abline(h = 0, lty = 2)
  if (filename$pulses == 1) {
    index = filename$tunits * (filename$time_pulse1 + time.delay)/filename$time_end
    arrows(x0 = filename$time[index], y0 = 0,
           x1 = filename$time[index], y1 = filename$charge[index],
           code = 3, length = 0.1, angle = 15)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$charge)/scale.factor - 0.10 * delta.y,
         substitute(paste("Q: ", q, " µC"), list(q = noquote(formatC(filename$charge[index], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
  } else {
    index1 = filename$tunits * (filename$time_pulse1 + time.delay)/filename$time_end
    index2 = filename$tunits * (filename$time_pulse2 + time.delay)/filename$time_end
    arrows(x0 = filename$time[index1], y0 = 0,
           x1 = filename$time[index1], y1 = filename$charge[index1],
           code = 3, length = 0.1, angle = 15)
    arrows(x0 = filename$time[index2], y0 = 0,
           x1 = filename$time[index2], y1 = filename$charge[index2],
           code = 3, length = 0.1, angle = 15)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$charge)/scale.factor - 0.10 * delta.y,
         substitute(paste(Q[f], ": ", qfor, " µC"), list(qfor = noquote(formatC(filename$charge[index1], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$charge)/scale.factor - 0.175 * delta.y,
         substitute(paste(Q[r], ": ", qrev, " µC"), list(qrev = noquote(formatC(filename$charge[index2], format = "f", digits = 2)))), 
         adj = c(0, NA), cex = 0.80)
    text(x = 0.8 * max(filename$time), 
         y = max(filename$charge)/scale.factor - 0.25 * delta.y,
         substitute(paste(Q[r]/Q[f], ": ", ratio), list(ratio = noquote(formatC(abs(filename$charge[index2]/filename$charge[index1]), format = "f", digits = 3)))), 
         adj = c(0, NA), cex = 0.80)
  }
}

#' diffusionGrid: function that displays the diffusion grid used in a cyclic voltammetry, a linear sweep voltammetry, a chronoamperometry, or a chronocoulometry simulation; increasing the scale.factor narrows the range of values on the x-axis so that diffusion profiles are displayed over a smaller range of distances

diffusionGrid = function(filename, 
                         species = c("Ox", "Red", "Z", "All"), 
                         scale.factor = 1){
  
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    stop("You need to install the plot3D package to use diffusionGrid.")
  }
  
  library(plot3D)
  
  species = match.arg(species)
  
#' in order to orient the diffusion grid so that rows correspond to times, with values on the y-axis and so that columns correspond to distances, with values on the x-axis, the transpose of the diffusion grid matrix is used  
  
  if (species == "Ox"){
    image2D(z = t(filename$oxdat), 
            y = filename$time, 
            x = filename$distance,
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Ox] (mM)", col = ramp.col(c("white", "blue")))
  } else if (species == "Red") {
    image2D(z = t(filename$reddata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Red] (mM)", col = ramp.col(c("white", "red")))
  } else if (species == "Z") {
    if (filename$mechanism == "E") {
     chem.conc = matrix (NA, nrow = filename$tunits, ncol = filename$xunits)
    } else {
      chem.conc = filename$chemdata
    }
    image2D(z = t(chem.conc), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)",
            ylim = c(0, filename$conc.bulk),
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Z] (mM)", col = ramp.col(c("white", "green")))
  } else if(species == "All") {
    old.par = par(mfrow = c(1, 3))
    image2D(z = t(filename$oxdata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)",
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Ox] (mM)", col = ramp.col(c("white", "blue")))
    image2D(z = t(filename$reddata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Red] (mM)", col = ramp.col(c("white", "red")))
    if (filename$mechanism == "E") {
      chem.conc = matrix (NA, nrow = filename$tunits, ncol = filename$xunits)
    } else {
      chem.conc = filename$chemdata
    }
    image2D(z = t(chem.conc), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Z] (mM)", col = ramp.col(c("white", "green")))
    par(old.par)
  }
}
