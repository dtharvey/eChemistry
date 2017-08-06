#' Notes: needs work on acccounting for chemical reactions; see cvTestProd.R for some initial work; need to look more closely at Gosser for general flow as he seems to do electrode (that is flux and current), then diffusion, then chemistry while work here is diffusion and then current with chemistry unclear; there seems to be a significant error in the paper by Brown for the EC mechansim in that the diffusion profile for the product is incorrect

#' To Do: break this file into parts for simulations, for basic plotting, for animations, and for other functions; continue exploring work in cvTestProd.R for ECr mechanism

#' The functions below simulate cyclic voltammetry (CV), linear-sweep voltammetry (LSV), both with and without stirring of the solution, and chronoamperometry (CA) experiments. Each simulation assumes an initial n-electron reduction or oxidation reaction and includes the option of an EC mechanism in which the product of the initial redox reaction undergoes an irreversible unimolecular chemical reaction. Noise, drawn at random from a normal distribution, is added to a simulation by specifiying the distribution's mean as 0 and its standard deviation as a percentage of the maximum current. The simulations use a set of diffusion grids to calculate and to store the concentrations of Ox and Red at discrete distances from the electrode surface at discrete times using the method outlined in Brown, J. H., J. Chem. Educ., 2015, 92, 1490-1496, and in Gosser, D. K. Cyclic Voltammetry Simulation and Analysis of Reaction Mechanisms, VCH, New York, 1993. Additional functions provide ways to display and to examine the results of these simulations.

#' The following is list of parameters used as inputs to the functions or as outputs of the functions

#' experiment: one of CV, LSV, CA, of CC
#' mechanism: one of E or EC
#' filetype: full (all data); reduced (subsample of data)
#' e.start: initial potential for all experiments (V)
#' e.switch: switching potential for cyclic voltammetry (V)
#' e.pulse: potential pulse for chronoamperometry (V)
#' e.end: final potential for linear-sweep voltammetry (V)
#' e.form: formal potential of the redox couple (V)
#' n: electrons for redox couple (unitless)
#' ko: standard heterogeneous electron transfer rate constant (cm/s)
#' kchem: homogeneous rate constant for chemical reaction (s^-1)
#' alpha: transfer coefficient (unitless)
#' d: diffusion coefficient for Ox and Red (cm^2/s)
#' area: surface area of the electrode (cm^2)
#' temp: temperature (K)
#' scan.rate: rate at which the potential is changed (V/s)
#' t.pulse1: time for first pulse for chronoamp (s)
#' t.pulse2: time for second pulse for chronoamp (s)
#' conc.bulk: concentration of analyte in bulk solution (mol/L)
#' oxdata: matrix giving [Ox] in diffusion grid (in mM)
#' reddata: matrix giving [Red] in diffusion grid (in mM)
#' t.units: number of time units for diffusion grid (unitless)
#' x.units: number of distance units for diffusion grid (unitless)
#' sd.noise: standard deviation as percent of imax (µA)
#' direction: -1 if Ox -> Red or +1 if Red -> Ox 

#' cycvolt: function to simulate a cyclic voltammetry experiment as either E only (kchem = 0) or an EC (kchem > 0); default conditions are from section 2.3.1 of Gosser's text, which is an E-only mechanism

cycvolt = function(e.start = 0.0, e.switch = -0.5, e.form = -0.25, 
                   n = 1, alpha = 0.50, ko = 1, kchem = 0,
                   d = 1e-5, area = 0.01, temp = 298.15, 
                   scan.rate = 1.0, conc.bulk = 1e-3, 
                   t.units = 1000, x.units = 100, sd.noise = 0) {
  
#' Test to ensure that t.units and x.units satisfy constraint that the number of distance units is less than (18 * number of time unit)^(0.5). The accuracy of the simulation improves with an increase in the number of time units and the number of distance units, but at the cost of an increase in the time needed to calculate the diffusion grid; the default condition of 1000 time units and 100 distance units has a total system time of approximately 0.1 s
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }
  
#' physical constants used in simulations: Faraday's contant (F) in C/mol and the gas constant (R) in J/K•mol
  
  f = 96485   
  r = 8.31451  
  
#' define the limits for the diffusion grid with respect to time and to distance, and calculatie additional simulation parameters
#'  t.tot: time to complete one full sweep from e.start to e.start (s)
#'  delta.t: increment in time (s)
#'  time: vector of discrete times for diffusion grid (s)
#'  x.tot: max distance chosen to exceed difusion limit (cm)
#'  delta.x: increment in distance (cm)
#'  distance: vector of discrete distances for diffusion grid (cm)
#'  lambda: a gathering of constants (unitless)
#'  direction: -1 for initial reduction and +1 for initial oxidation
#'  cox.bulk: bulk concentration of Ox (mol/cm^3)
#'  cred.bulk: bulk concentration of Red (mol/cm^3)
  
  t.tot = 2 * abs((e.start - e.switch))/scan.rate
  delta.t = t.tot/t.units
  time = seq(0, t.tot, delta.t)
  x.tot = 6 * sqrt(d * t.tot)
  delta.x = x.tot/x.units
  distance = seq(0, x.tot, delta.x)
  lambda = d * delta.t/(delta.x)^2
  if (e.start > e.switch) {
    direction = -1
    cox.bulk = conc.bulk
    cred.bulk = 0
  } else {
    direction = +1
    cox.bulk = 0
    cred.bulk = conc.bulk
  }
  
#' create vector of discrete applied potentials
  
  pot_forward = seq(e.start, e.switch, direction * scan.rate * delta.t)
  pot_reverse = seq(e.switch, e.start, -direction * scan.rate * delta.t)
  potential = c(pot_forward, pot_reverse[-1])
  
# calculatie the potential-dependent forward (kf) and reverse (kb) heterogeneous electron-transfer rate constants

  kf = ko * exp(-alpha * n * f * (potential - e.form)/(r*temp))
  kb = ko * exp((1 - alpha) * n * f * (potential - e.form)/(r*temp))
  
#' initialize diffusion grid (rows = time; cols = distance) using bulk concentrations for Ox and for Red and adjusting concentrations to mol/cm^3; the actual concentrations are calculated later
  
  dif.ox = matrix(cox.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.red = matrix(cred.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
#' create vectors for fluxes and current, which are calculated later; the initial values here are not important as actual values are calculated later
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
#' calculate diffusion grid over time and, at each time, over distance; for each time the diffusion grid for Ox and for Red first is calculated at all distances except for that at the electrode surface; next, for each time, the flux of Ox and Red to the electrode surface are used to calculate their concentrations at the electrode surface; and, finally, for each time, the current is calculated
  
  for (i in 2:(t.units + 1)){
    for (j in 2:x.units) {
      if (direction == -1) {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kchem * delta.t * dif.red[i-1, j]
      } else {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kchem * delta.t * dif.ox[i-1, j]
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
      }
    }
    jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
    jred[i] = -jox[i]
    dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
    dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
    current.total[i] = -n * f * area * jox[i]
  }
  
#' add noise to the current; note default is sd.noise = 0, which returns the pure, noise-free simulated cyclic voltammogram

  noise = rnorm(t.units + 1, mean = 0, 
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
#' return original inputs and calculated results as a list for use with other functions

  if (kchem > 0) {mechanism = "EC"} else {mechanism = "E"}
  
  output = list("expt" = "CV",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox*10^6, 
                "reddata" = dif.red*10^6, 
                "formalE" = e.form,
                "initialE" = e.start,
                "switchE" = e.switch,
                "electrons" = n,
                "k_o" = ko,
                "kchem" = kchem,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "scanrate" = scan.rate,
                "conc_bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction
  )

  invisible(output)
}

#' linsweep: function to simulate linear sweep voltammetry from an initial potential to an end potential; parameters are identical to those in the function cycvolt with the expection that e.switch is replaced with e.end and with addition of stir.rate, which takes values of off, slow, medium, or fast; see comments for the function cycvolt for additional details not noted here

linsweep = function(e.start = 0.0, e.end = -1, e.form = -0.25, 
                    n = 1, alpha = 0.50, ko = 1, kchem = 0,
                    d = 1e-5, area = 0.01, temp = 298.15, 
                    scan.rate = 1.0, conc.bulk = 1e-3, 
                    t.units = 1000, x.units = 100, 
                    stir.rate = "off", sd.noise = 0) {
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }

#' translate stir.rate to a factor that limits the thickness of the diffusion layer when stirring is turned on; this is accomplished in the loops for calculating the diffusion grid by dividing x.units by the scaling.factor defined here 

  if (stir.rate == "off") {
    scaling.factor = 1
  } else if (stir.rate == "slow") {
    scaling.factor = 10
  } else if (stir.rate == "medium") {
    scaling.factor = 25
  } else if (stir.rate == "fast") {
    scaling.factor = 50
  } else {
    stop("Options for stir rate are off, slow, medium, or fast.")
  }
  
  f = 96485  
  r = 8.31451  
  
  t.tot = abs((e.start - e.end))/scan.rate
  delta.t = t.tot/t.units
  time = seq(0, t.tot, delta.t)
  x.tot = 6 * sqrt(d * t.tot)
  delta.x = x.tot/x.units
  distance = seq(0, x.tot, delta.x)
  lambda = d * delta.t/(delta.x)^2
  if (e.start > e.end) {
    direction = -1
    cox.bulk = conc.bulk
    cred.bulk = 0
  } else {
    direction = +1
    cox.bulk = 0
    cred.bulk = conc.bulk
  }
  
  potential = seq(e.start, e.end, direction * scan.rate * delta.t) 
  
  kf = ko * exp(-alpha * n * f * (potential - e.form)/(r*temp))
  kb = ko * exp((1 - alpha) * n * f * (potential - e.form)/(r*temp))
  
  dif.ox = matrix(cox.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.red = matrix(cred.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
  for (i in 2:(t.units + 1)){
    for (j in 2:(x.units/scaling.factor)) {
      if (direction == -1) {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kchem * delta.t * dif.red[i-1, j]
      } else {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kchem * delta.t * dif.ox[i-1, j]
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
      }
    }
    jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
    jred[i] = -jox[i]
    dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
    dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
    current.total[i] = -n * f * area * jox[i]
  }
  
  noise = rnorm(t.units + 1, mean = 0,
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
  if (kchem > 0) {mechanism = "EC"} else {mechanism = "E"}
  
  output = list("expt" = "LSV",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox*10^6, 
                "reddata" = dif.red*10^6, 
                "formalE" = e.form,
                "initialE" = e.start,
                "endE" = e.end,
                "electrons" = n,
                "k_o" = ko,
                "kchem" = kchem,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "scanrate" = scan.rate,
                "conc_bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction,
                "stir_rate" = stir.rate
  )
  invisible(output)
}

#' chronoamp: function to simulate chronoamperometry; parameters are similar to those in the function cycvolt with the addition of pulses, which takes a value of 1 for a single pulse experiment and a value of 2 for a double pulse experiment; the total elapsed time to the first pulse is given by t.1, the total elapsed time to the second pulse, if used, is t.2, and the total elapsed time to the experiment's end is t.end; see comments for the function cycvolt for additional details not noted here

chronoamp = function(e.start = 0.0, e.pulse = -0.5, e.form = -0.25,
                     pulses = 1, t.1 = 10, t.2 = 0, t.end = 30,
                     n = 1, alpha = 0.50, ko = 1, kchem = 0,
                     d = 1e-5, area = 0.01, temp = 298.15, 
                     conc.bulk = 1e-3, 
                     t.units = 1000, x.units = 100, sd.noise = 0) {
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }
  
  f = 96485 
  r = 8.31451 
  
  t.tot = t.end
  delta.t = t.tot/t.units
  time = seq(0, t.tot, delta.t)
  x.tot = 6 * sqrt(d * t.tot)
  delta.x = x.tot/x.units
  distance = seq(0, x.tot, delta.x)
  lambda = d * delta.t/(delta.x)^2
  if (e.start > e.pulse) {
    direction = -1
    cox.bulk = conc.bulk
    cred.bulk = 0
  } else {
    direction = +1
    cox.bulk = 0
    cred.bulk = conc.bulk
  }
  
  if (pulses == 1) { 
    pot1 = rep(e.start, round(t.units * t.1/t.tot, digits = 0))
    pot2 = rep(e.pulse, t.units - length(pot1) + 1)
    potential = c(pot1, pot2)
  } else {
    pot1 = rep(e.start, round(t.units * t.1/t.tot, digits = 0))
    pot2 = rep(e.pulse, round(t.units * (t.2 - t.1)/t.tot, digits = 0))
    pot3 = rep(e.start, t.units - length(pot1) - length(pot2) + 1)
    potential = c(pot1, pot2, pot3)
  }
  
  kf = ko * exp(-alpha * n * f * (potential - e.form)/(r*temp))
  kb = ko * exp((1 - alpha) * n * f * (potential - e.form)/(r*temp))
  
  dif.ox = matrix(cox.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.red = matrix(cred.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
  for (i in 2:(t.units + 1)){
    for (j in 2:x.units) {
      if (direction == -1) {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kchem * delta.t * dif.red[i-1, j]
      } else {
        dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kchem * delta.t * dif.ox[i-1, j]
        dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
      }
    }
    jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
    jred[i] = -jox[i]
    dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
    dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
    current.total[i] = -n * f * area * jox[i]
  }
  
  noise = rnorm(t.units + 1, mean = 0,
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
  if (kchem > 0) {mechanism = "EC"} else {mechanism = "E"}
  
  output = list("expt" = "CA",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox*10^6, 
                "reddata" = dif.red*10^6, 
                "formalE" = e.form,
                "initialE" = e.start,
                "pulseE" = e.pulse,
                "electrons" = n,
                "k_o" = ko,
                "kchem" = kchem,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "conc_bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction,
                "pulses" = pulses,
                "time_pulse1" = t.1,
                "time_pulse2" = t.2,
                "time_end" = t.end
  )
  invisible(output)
}

#' chronocoul: function that converts the result of a chronoamperometry simulation into its corresponding chonocoulometry simulation by integrating current over time using trapazoidal rule integration (adapted here from the function trapz in the caTools package)

chronocoul = function(filename){
  
  if (filename$expt != "CA") {
    stop("This file is not from a chronoamperometry simulation.")
  }
  
  charge = rep(0, length(filename$current))
  
  for (i in 2:length(charge)) {
    x = 2:i
    charge[i] = as.double((filename$time[x] - filename$time[x-1]) %*% (filename$current[x] + filename$current[x - 1])/2)
  }
  
  output = list("expt" = "CC",
                "mechanism" = filename$mechanism,
                "file_type" = "full",
                "charge" = charge, 
                "potential" = filename$potential,
                "time" = filename$time, 
                "distance" = filename$distance, 
                "oxdata" = filename$oxdata, 
                "reddata" = filename$reddata, 
                "formalE" = filename$formalE,
                "initialE" = filename$initialE,
                "pulseE" = filename$pulseE,
                "electrons" = filename$electrons,
                "k_o" = filename$k_o,
                "kchem" = filename$kchem,
                "alpha" = filename$alpha,
                "diffcoef" = filename$diffcoef,
                "area" = filename$area,
                "temperature" = filename$temperature,
                "conc_bulk" = filename$conc_bulk,
                "tunits" = filename$tunits,
                "xunits" = filename$xunits,
                "sdnoise" = filename$sdnoise,
                "direction" = filename$direction,
                "pulses" = filename$pulses,
                "time_pulse1" = filename$time_pulse1,
                "time_pulse2" = filename$time_pulse2,
                "time_end" = filename$time_end
  )
  
  invisible(output)

}

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
  
  ymax = filename$conc_bulk * 1000
  
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
    legend(x = "topright", legend = c("OX", "RED"), 
           fill = c("blue", "red"), 
           bty = "n", inset = c(0.05, 0.05))
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
          y = filename$proddata[index, ], 
          lwd = 3, col = "green")
    legend(x = "topright", legend = c("OX", "RED", "PROD"), 
           fill = c("blue", "red", "green"), 
           bty = "n", inset = c(0.05, 0.05, 0.05))
    grid()
  } 
}

#' plotGrid: function that plots eight diffusion profiles---at times that are 10%, 20%, 30%, 40%, 60%, 70%, 80%, and 90% of the total time for the experiment---around a central plot that shows the voltammogram or chronoamperogram from an object created useig cycvolt, linsweep, or chronoamp; all plots have default titles

plotGrid = function(filename) {
  
  #' adjust the graphical parameters to create a 3 by 3 grid of graphical windows and to adjust the margins for each graphical window
  
  par(mfrow = c(3,3), mar = c(4.1, 4.1, 2.1, 2.1))
  
  #' set the times for the diffusion profiles 
  
  t = round(c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9) * (filename$tunits), digits = 0)
  
  
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

#' plotVgram: function that plots 1--5 voltammograms (cyclic or linear sweep) on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; line widths, line types, and line colors have default values that can be adjusted

plotVgram = function(filenames = list(file1, file2, file3, file, file5), 
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
    plot(x = filenames[[1]]$potential[20], 
         y = filenames[[1]]$current[1],
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

#' plotCAgram: function that plots 1--5 chronoamperograms on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; setting scale to a value less than 1 adjusts the y-axis limits so that the limits are not set by the current spikes; line widths, line types, and line colors have default values that can be adjusted

plotCAgram = function(filenames = list(file1, file2, file3, file, file5),
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
    plot(x = filenames[[1]]$time[20], 
         y = filenames[[1]]$current[1],
         xlim = c( xmin[xmin_id], xmax[xmax_id]), 
         ylim = c(ymin[ymin_id], ymax[ymax_id]),
         type = "p", lwd = line_widths[1] , lty = line_types[1], 
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

#' plotCCgram: function that plots 1--5 chronocoulomograms on a single set of axes; files are passed as a list; the default plot does not include a legend or a title, but providing a vector of character strings to legend_text adds a legend to the final plot and adding a character string for main_title adds a title to the plot; setting scale to a value less than 1 adjusts the y-axis limits so that the limits are not set by the current spikes; line widths, line types, and line colors have default values that can be adjusted

plotCCgram = function(filenames = list(file1, file2, file3, file, file5),
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
    plot(x = filenames[[1]]$time[20], 
         y = filenames[[1]]$charge[1],
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

#' annotateCV: function that plots a cyclic voltammogram and annotates it with values for Epc, Epa, delta E, ip,c, ip,a, and the peak current ratio; the percentage of points to include in determining the peak currents is set using cathodic.per and anodic.per; to be measurable, the return peak current must satisfy a minimum threshold value and must be greater than the predicted baseline current 

annotateCV = function(filename, forward.per = 5, reverse.per = 5, threshold = 0.05) {
  
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
  
#' create linear models to predict baseline currents using parameters forward.per and reverse.per; lm1 models the baseline for the foward reaction and uses the first X% of the data points to model the baseline and a first-order relationship between current and potential; lm2 models the baseline for the reverse reaction and uses the first X% of the data points beginning with the switching potential to predict the current decay in the absence of a change in potential using a model in which the current decays as function of t^(-0.5)
  
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
  
#' annotations for anodic peak potential and current
  
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
  
#' annotations for delta E, Eavg, and current ratio
  
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

#' annotateLSV: function that plots a linear-sweep voltammogram and annotates it with values for either Epc and ip,c, or for Epa and ip,a,; the percentage of points to include in determining the peak current is set using potential.per; see annotateCV for comments on the code

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

#' plotlyVgram: function that plots an interactive cyclic voltammogram or linear sweep voltammogram using the plotly package; when hovering the cursor near the data it displays values for the current and the potential

plotlyVgram = function(filename) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("You need to install the plotly package to use plotlyVgram.")
  }
  library(plotly)
  
  plot_ly(x = filename$potential, y = filename$current, type = "scatter", mode = "lines") %>% layout(xaxis = list(autorange = "reversed", title = "potential (V)", showspikes = TRUE)) %>% layout(yaxis = list(title = "current (µA)", showspikes = TRUE)) %>% layout(hovermode = 'closest')
}

#' plotlyCAgram: function that plots an interactive chronoamperogram using the plotly package; when hovering the cursor near the data it displays values for the current and the potential

plotlyCAgram = function(filename, scale.factor = 1) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("You need to install the plotly package to use plotlyCAgram.")
  }
  library(plotly)
  
  plot_ly(x = filename$time, y = filename$current, type = "scatter", mode = "lines") %>% layout(xaxis = list(title = "time (s)", showspikes = TRUE)) %>% layout(yaxis = list(title = "current (µA)", showspikes = TRUE)) %>% layout(hovermode = 'closest')
}

#' plotlyCCgram: function that plots an interactive chronocoulogram using the plotly package; when hovering the cursor near the data it displays values for the current and the potential

plotlyCCgram = function(filename, scale.factor = 1) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("You need to install the plotly package to use plotlyCAgram.")
  }
  library(plotly)
  
  plot_ly(x = filename$time, y = filename$charge, type = "scatter", mode = "lines") %>% layout(xaxis = list(title = "time (s)", showspikes = TRUE)) %>% layout(yaxis = list(title = "charge (µC)", showspikes = TRUE)) %>% layout(hovermode = 'closest')
}

#' tableVgram: function that produces a table of times, potentials, currents, and concentrations of Ox and of Red at the electrode surface for a cyclic voltammogram or a linear sweep voltammogram; requires the DT and the magrittr packages

tableVgram = function(filename) {
  
  #' check to see if DT and magritrr packages are installed
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("You need to install the DT package to use tableVgram.")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("You need to install the magrittr package to use tableVgram.")
  }
  
  library(DT)
  library(magrittr)
  
#' create data frame of values and assign column names
  df = data.frame(filename$time, filename$potential, filename$current, 
                  filename$oxdata[ , 1], filename$reddata[ , 1])
  colnames(df) = c("t (s)", "E (V)", "i (µA)", 
                   "[Ox]x=0 (mM)", "[Red]x=0 (mM)")
  
# use the DT datatable command to create HTML table
  datatable(df, options = list(columnDefs = list(list(className = 'dt-center', targets = c(1:5))))) %>% 
    formatRound("t (s)", 3) %>% 
    formatRound("E (V)", 3) %>%
    formatRound("i (µA)", 3) %>%
    formatRound("[Ox]x=0 (mM)", 3) %>%
    formatRound("[Red]x=0 (mM)", 3)
}

#' tableCAgram: function that produces a table of times, potentials, currents, and concentrations of Ox and of Red at the electrode surface for a chronoamperogram; requires the DT and the magrittr packages

tableCAgram = function(filename) {
  
#' check to see if DT and magritrr packages are installed
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("You need to install the DT package to use tableVgram.")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("You need to install the magrittr package to use tableVgram.")
  }
  
  library(DT)
  library(magrittr)
  
#' create data frame of values and assign column names
  df = data.frame(filename$time, filename$potential, filename$current, 
                  filename$oxdata[ , 1], filename$reddata[ , 1])
  colnames(df) = c("t (s)", "E (V)", "i (µA)", 
                   "[Ox]x=0 (mM)", "[Red]x=0 (mM)")
  
# use the DT datatable command to create HTML table
  datatable(df, options = list(columnDefs = list(list(className = 'dt-center', targets = c(1:5))))) %>% 
    formatRound("t (s)", 3) %>% 
    formatRound("E (V)", 3) %>%
    formatRound("i (µA)", 3) %>%
    formatRound("[Ox]x=0 (mM)", 3) %>%
    formatRound("[Red]x=0 (mM)", 3)
}

#' tableCCgram: function that produces a table of times, potentials, charges, and concentrations of Ox and of Red at the electrode surface for a chronocoulogram; requires the DT and the magrittr packages

tableCCgram = function(filename) {
  
  #' check to see if DT and magritrr packages are installed
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("You need to install the DT package to use tableVgram.")
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop("You need to install the magrittr package to use tableVgram.")
  }
  
  library(DT)
  library(magrittr)
  
#' create data frame of values and assign column names
  df = data.frame(filename$time, filename$potential, filename$charge, 
                  filename$oxdata[ , 1], filename$reddata[ , 1])
  colnames(df) = c("t (s)", "E (V)", "Q (µC)", 
                   "[Ox]x=0 (mM)", "[Red]x=0 (mM)")
  
# use the DT datatable command to create HTML table
  datatable(df, options = list(columnDefs = list(list(className = 'dt-center', targets = c(1:5))))) %>% 
    formatRound("t (s)", 3) %>% 
    formatRound("E (V)", 3) %>%
    formatRound("Q (µC)", 3) %>%
    formatRound("[Ox]x=0 (mM)", 3) %>%
    formatRound("[Red]x=0 (mM)", 3)
}

#' sampleVgram: used to create a reduced data file of potentials and currents for a cyclic voltammogram or a linear sweep voltammogram; data.reduction gives the percentage of data to keep, which then is evenly spaced using a calculated increment 

sampleVgram = function(filename, data.reduction = 1) {
  potential = filename$potential
  current = filename$current
  len = length(filename$potential)
  delta.data = len * data.reduction/100
  pot.new = potential[seq(1, len, len/delta.data)]
  cur.new = current[seq(1, len, len/delta.data)]
  if (filename$expt == "CV") {
    output = list("expt" = filename$expt,
                  "mechanism" = filename$mechanism,
                  "file_type" = "reduced",
                  "current" = cur.new, 
                  "potential" = pot.new,
                  "formalE" = filename$formalE,
                  "initialE" = filename$initialE,
                  "switchE" = filename$switchE,
                  "electrons" = filename$electrons,
                  "k_o" = filename$k_o,
                  "kchem" = filename$kchem,
                  "alpha" = filename$alpha,
                  "diffcoef" = filename$diffcoef,
                  "area" = filename$area,
                  "temperature" = filename$temperature,
                  "scanrate" = filename$scanrate,
                  "conc_bulk" = filename$conc_bulk,
                  "tunits" = filename$tunits,
                  "xunits" = filename$xunits,
                  "sdnoise" = filename$sdnoise,
                  "direction" = filename$direction,
                  "increment" = delta.data
    )
    invisible(output)
  } else {
    output = list("expt" = filename$expt,
                  "mechanism" = filename$mechanism,
                  "file_type" = "reduced",
                  "current" = cur.new, 
                  "potential" = pot.new,
                  "formalE" = filename$formalE,
                  "initialE" = filename$initialE,
                  "endE" = filename$endE,
                  "electrons" = filename$electrons,
                  "k_o" = filename$k_o,
                  "kchem" = filename$kchem,
                  "alpha" = filename$alpha,
                  "diffcoef" = filename$diffcoef,
                  "area" = filename$area,
                  "temperature" = filename$temperature,
                  "scanrate" = filename$scanrate,
                  "conc_bulk" = filename$conc_bulk,
                  "tunits" = filename$tunits,
                  "xunits" = filename$xunits,
                  "sdnoise" = filename$sdnoise,
                  "direction" = filename$direction,
                  "stir_rate" = filename$stir_rate,
                  "increment" = data.reduction
    )
    invisible(output)
  }
}

#' sampleCAgram: used to create a reduced data file of times and currents for a chronoamperogram; data.reduction gives the percentage of data to keep, which then is evenly spaced using a calculated increment 

sampleCAgram = function(filename, data.reduction = 1) {
  time = filename$time
  current = filename$current
  len = length(filename$time)
  delta.data = len * data.reduction/100
  time.new = time[seq(1, len, len/delta.data)]
  cur.new = current[seq(1, len, len/delta.data)]
  output = list("expt" = filename$expt,
                "mechanism" = filename$mechanism,
                "file_type" = "reduced",
                "current" = cur.new, 
                "time" = time.new,
                "formalE" = filename$formalE,
                "initialE" = filename$initialE,
                "pulseE" = filename$pulseE,
                "electrons" = filename$electrons,
                "k_o" = filename$k_o,
                "kchem" = filename$kchem,
                "alpha" = filename$alpha,
                "diffcoef" = filename$diffcoef,
                "area" = filename$area,
                "temperature" = filename$temperature,
                "conc_bulk" = filename$conc_bulk,
                "tunits" = filename$tunits,
                "xunits" = filename$xunits,
                "sdnoise" = filename$sdnoise,
                "direction" = filename$direction,
                "increment" = delta.data,
                "pulses" = filename$pulses,
                "time_pulse1" = filename$time_pulse1,
                "time_pulse2" = filename$time_pulse2,
                "time_end" = filename$time_end
    )
    invisible(output)
}

#' vgramGIF: function that creates an animated gif that consists of 41 individual plots of a cyclic voltammogram or linear sweep voltammogram

vgramGIF = function(filename, gif.name = "vgram"){
  
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

#' cagramGIF: function that creates an animated gif that consists of 41 individual plots of a chronoamperogram

cagramGIF = function(filename, gif.name = "cagram"){
  
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

#' ccgramGIF: function that creates an animated gif that consists of 41 individual plots of a chronocoulogram

ccgramGIF = function(filename, gif.name = "ccgram"){
  
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

#' diffusionGIF: function that creates an animated gif that consists of 41 individual plots of a diffusion profile for a cyclic voltammetry, linear sweep voltammetry, or chronoamperometry experiment

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
      }
    }}, movie.name = paste0(gif.name,".gif"))
  
#' reset the original values for ani.options   
  ani.options(old.ani)
}

#' cvGIF: function that creates an animated gif that consists of 41 individual plots of a cyclic voltammetry experiment, including both the cyclic voltammogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' lsvGIF: function that creates an animated gif that consists of 41 individual plots of a linear sweep voltammetry experiment, including both the linear sweep voltammogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

# caGIF: function that creates an animated gif that consists of 41 individual plots of a chronoamperometry experiment, including both the chronoamperogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

# ccGIF: function that creates an animated gif that consists of 41 individual plots of a chronocoulometry experiment, including both the chronocoulogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' vgramHTML: function that creates an HTML file with a movie of a cyclic voltammogram or a linear sweep voltammogram

vgramHTML = function(filename, html.name = "vgram"){
  
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

#' cagramHTML: function that creates an HTML file with a movie of a chronoamperogram

cagramHTML = function(filename, html.name = "cagram"){
  
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

#' ccgramHTML: function that creates an HTML file with a movie of a chronocoulogram

ccgramHTML = function(filename, html.name = "ccgram"){
  
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

#' diffusionHTML: function that creates an HTML file with a movie of a diffusion profile a diffusion profile for a cyclic voltammetry, linear sweep voltammetry, or chronoamperometry experiment

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
           ylim = c(0, 1050 * filename$conc_bulk),
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        lines(x = filename$distance, y = filename$reddata[i, ], 
              lwd = 3, col = "red")
        legend(x = "right", legend = c("Ox", "Red"), 
               fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.01))
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

#' cvHTML: function that creates an animated HTML file that consists of 41 individual plots of a cyclic voltammetry experiment, including both the cyclic voltammogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' lsvHTML: function that creates an animated HTML file that consists of 41 individual plots of a linear sweep voltammetry experiment, including both the linear sweep voltammogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' caHTML: function that creates an animated HTML that consists of 41 individual plots of a chronoamperometry experiment, including both the chronoamperogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' ccHTML: function that creates an animated HTML that consists of 41 individual plots of a chronocoulometry experiment, including both the chronocoulogram and the corresponding diffuions profiles

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
           ylim = c(0, 1050 * filename$conc_bulk), 
           xlab = "distance from electrode (cm)", 
           ylab = "concentration (mM)")
      grid()
      lines(x = filename$distance, y = filename$reddata[i, ], 
            lwd = 3, col = "red")
      if (filename$mechanism == "EC") {
        lines(x = filename$distance, 
              y = 1000 * filename$conc_bulk - filename$oxdata[i, ] - filename$reddata[i, ], 
              lwd = 3, col = "green")
        legend(x = "right", legend = c("OX", "RED", "PROD"), 
               fill = c("blue", "red", "green"), 
               bty = "n", inset = c(0.05, 0.05, 0.05))
      } else {
        legend(x = "right", legend = c("OX", "RED"), fill = c("blue", "red"),
               bty = "n", inset = c(0.05, 0.05))
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

#' diffusionGrid: function that displays the diffusion grid used in a cyclic voltammetry, a linear sweep voltammetry, a chronoamperometry, or a chronocoulometry simulation; increasing the scale.factor narrows the range of values on the x-axis so that diffusion profiles are displayed over a smaller range of distances

diffusionGrid = function(filename, species = "Ox", scale.factor = 1){
  
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    stop("You need to install the plot3D package to use diffusionGrid.")
  }
  library(plot3D)
  
#' in order to orient the diffusion grid so that rows correspond to times, with values on the y-axis and so that columns correspond to distances, with values on the x-axis, the transpose of the diffusion grid matrix is used  
  
  if (species == "Ox"){
    image2D(z = t(filename$oxdat), 
            y = filename$time, 
            x = filename$distance,
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Ox] (mM)", col = ramp.col(c("white", "darkblue")))
  } else if (species == "Red") {
    image2D(z = t(filename$reddata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Red] (mM)", col = ramp.col(c("white", "darkblue")))
  } else if(species == "both") {
    old.par = par(mfrow = c(1, 2))
    image2D(z = t(filename$oxdata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)",
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Ox] (mM)", col = ramp.col(c("white", "darkblue")))
    image2D(z = t(filename$reddata), 
            y = filename$time, 
            x = filename$distance, 
            contour = FALSE, 
            ylab = "time (s)", xlab = "distance (cm)", 
            xlim = c(min(filename$distance/scale.factor), 
                     max(filename$distance/scale.factor)),
            clab = "[Red] (mM)", col = ramp.col(c("white", "darkblue")))
    par(old.par)
  }
}

randsev = function(filename) {
  f = 96485   
  r = 8.31451
  temp = filename$temperature
  diffcoef = filename$diffcoef
  area = filename$area
  n = filename$electrons
  conc = filename$conc_bulk
  scan.rate = filename$scanrate
  pred.current = 0.446 * f * area * n * conc * sqrt((n * f * diffcoef *scan.rate)/(r * temp)) * 1000
  if (filename$direction == -1) {
    meas.current = max(filename$current)
  } else {
    meas.current = -min(filename$current)
  }
  out = list("measured_current" = meas.current, "predicted_current" = pred.current)
}
