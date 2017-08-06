#' The functions below simulate cyclic voltammetry (cvSim), linear-sweep voltammetry (lsvSim), both with and without stirring of the solution, chronoamperometry (caSim), and chronocoulometry (ccSim) experiments. Each simulation assumes an n-electron rdeox reaction and includes options for an EC mechanism, in which the product of the initial redox reaction undergoes a unimolecular chemical reaction, or a CE mechanism, in which an initial unimolecular chemical reaction generates the electroactive species. Noise from a random normal distribution is included by specifiying the distribution's mean as 0 and its standard deviation as a percentage of the maximum current. Each simulations uses a set of diffusion grids to calculate and to store the concentrations of each species at discrete distances from the electrode surface at discrete times using the method outlined in Brown, J. H., J. Chem. Educ., 2015, 92, 1490-1496, and in Gosser, D. K. Cyclic Voltammetry Simulation and Analysis of Reaction Mechanisms, VCH, New York, 1993. Additional functions provide ways to display and to examine the results of these simulations.

#' The following is a list of parameters used as inputs to the functions or as outputs of the functions

#' experiment: one of CV, LSV, CA, of CC
#' mechanism: one of E, EC, or CE
#' filetype: full (all data); reduced (subsample of data)
#' e.start: initial potential for all experiments (V)
#' e.switch: switching potential for cyclic voltammetry (V)
#' e.pulse: potential pulse for chronoamperometry (V)
#' e.end: final potential for linear-sweep voltammetry (V)
#' e.form: formal potential of the redox couple (V)
#' n: electrons for redox couple (unitless)
#' ko: standard heterogeneous electron transfer rate constant (cm/s)
#' kcf, kcr: homogeneous rate constants for chemical reaction (s^-1)
#' alpha: transfer coefficient (unitless)
#' d: diffusion coefficient for Ox and Red (cm^2/s)
#' area: surface area of the electrode (cm^2)
#' temp: temperature (K)
#' scan.rate: rate at which the potential is changed (V/s)
#' t.pulse1: time for first pulse for chronoamp (s)
#' t.pulse2: time for second pulse for chronoamp (s)
#' conc.bulk: total concentration of Ox, Red, and Chem (mol/L)
#' oxdata: matrix giving [Ox] in diffusion grid (in mM)
#' reddata: matrix giving [Red] in diffusion grid (in mM)
#' chem.data: matrix giving [Chem] in diffusion grid (in mM)
#' t.units: number of time units for diffusion grid (unitless)
#' x.units: number of distance units for diffusion grid (unitless)
#' sd.noise: standard deviation as percent of imax (µA)
#' direction: -1 if Ox -> Red or +1 if Red -> Ox 

#' cvSim: function to simulate a cyclic voltammetry experiment as either an E, EC, or CE mechanism, where E is a redox reaction and where C is a chemical reaction that either precedes or follows the redox reaction; the default conditions are from section 2.3.1 of Gosser's text, which is an E mechanism

cvSim = function(e.start = 0.0, e.switch = -0.5, e.form = -0.25,
                 mechanism = c("E", "EC", "CE"),
                 ko = 1, kcf = 0, kcr = 0,
                 n = 1, alpha = 0.50, d = 1e-5, area = 0.01, 
                 temp = 298.15, scan.rate = 1.0, conc.bulk = 1e-3,
                 t.units = 1000, x.units = 100, sd.noise = 0) {
  
#' Test to ensure that t.units and x.units satisfy constraint that the number of distance units is less than (18 * number of time unit)^(0.5). The accuracy of the simulation improves with an increase in the number of time units and the number of distance units, but at the cost of an increase in the time needed to calculate the diffusion grid; the default condition of 1000 time units and 100 distance units has a total system time of approximately 0.1 s
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }
  
  mechanism = match.arg(mechanism)
  
  if (mechanism == "E") {
    if (kcf != 0 | kcr != 0) {
      stop("For the E mechanism, kcf and kcr must have values of 0.")
    }
  }
  
  if (mechanism == "CE") {
    if (kcf <= 0) {
      stop("For the CE mechanism, kcf must have a value greater than 0.")
    }
  }
  
#' physical constants used in simulations: Faraday's contant (F) in C/mol and the gas constant (R) in J/K•mol
  
  f = 96485   
  r = 8.31451  
  
#' define the limits for the diffusion grid with respect to time and to distance, and calculate additional simulation parameters, including bulk concentrations of all species
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
  
  t.tot = 2*abs((e.start - e.switch))/scan.rate
  delta.t = t.tot/t.units
  time = seq(0, t.tot, delta.t)
  x.tot = 6 * sqrt(d * t.tot)
  delta.x = x.tot/x.units
  distance = seq(0, x.tot, delta.x)
  lambda = d * delta.t/(delta.x)^2
  if (e.start > e.switch) {
    direction = -1
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = conc.bulk - cchem.bulk
      cred.bulk = 0
    } else {
      cox.bulk = conc.bulk
      cred.bulk = 0
      cchem.bulk = 0
    }
  } else {
    direction = +1
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = 0
      cred.bulk = conc.bulk - cchem.bulk
    } else {
      cox.bulk = 0
      cred.bulk = conc.bulk
      cchem.bulk = 0
    }
  }
  
#' check to ensure that the number of time units satisfies Gosser's requirement for accuracy when using an EC or CE mechanism
  
  if (mechanism != "E") {
    min_tunits = 4 * t.tot * kcf
    if (t.units < min_tunits) {
      stop(paste("Minimum time units is", min_tunits, "if kcf =", kcf, "and with a total scan time of", t.tot, "s."))
    }
  }
  
#' create vector of discrete applied potentials
  
  pot_forward = seq(e.start, e.switch, direction * scan.rate * delta.t)
  pot_reverse = seq(e.switch, e.start, -direction * scan.rate * delta.t)
  potential = c(pot_forward, pot_reverse[-1])
  
# calculate the potential-dependent forward (kf) and reverse (kb) heterogeneous electron-transfer rate constants
  
  kf = ko * exp(-alpha * n * f * (potential - e.form)/(r*temp))
  kb = ko * exp((1 - alpha) * n * f * (potential - e.form)/(r*temp))
  
#' initialize diffusion grids (rows = time; cols = distance) using bulk concentrations for Ox, Red, and Chem and adjusting concentrations to mol/cm^3; the actual concentrations are calculated later
  
  dif.ox = matrix(cox.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.red = matrix(cred.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.chem = matrix(cchem.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
#' create vectors for fluxes and current, which are calculated later; the initial values here are not important as actual values are calculated later
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
#' calculate diffusion grids over time and, at each time, over distance; for each time the diffusion grids is calculated at all distances except for at the electrode surface; next, for each time, the flux of each species to the electrode surface is used to calculate their concentrations at the electrode surface; and, finally, for each time, the current is calculated
  
  if (mechanism == "CE") {
    for (i in 2:(t.units + 1)){
      for (j in 2:x.units) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) + kcf * delta.t * dif.chem[i-1, j] - kcr * delta.t * dif.ox[i-1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.ox[i-1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) 
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) + kcf * delta.t*dif.chem[i-1, j] - kcr * delta.t * dif.red[i-1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.red[i-1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  } else {
    for (i in 2:(t.units + 1)){
      for (j in 2:x.units) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kcf * delta.t * dif.red[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.red[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kcf * delta.t * dif.ox[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.ox[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  }
  
#' if desired, add noise to the current; note default is sd.noise = 0, which returns the pure, noise-free simulated cyclic voltammogram
  
  noise = rnorm(t.units + 1, mean = 0, 
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
#' return original inputs and calculated results as a list for use with other functions
  
  output = list("expt" = "CV",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox * 10^6, 
                "reddata" = dif.red * 10^6, 
                "chemdata" = dif.chem * 10^6,
                "formalE" = e.form,
                "initialE" = e.start,
                "switchE" = e.switch,
                "electrons" = n,
                "ko" = ko,
                "kcf" = kcf,
                "kcr" = kcr,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "scanrate" = scan.rate,
                "conc.bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction
  )
  
  invisible(output)
}

#' lsvSim: function to simulate a linear sweep voltammetry experiment from an initial potential to an end potential; parameters are identical to those in the function cvSim with the expection that e.switch is replaced with e.end and with addition of stir.rate, which takes values between 0 and 1; see comments for the function cvSim for additional details not noted here

lsvSim = function(e.start = 0.0, e.end = -1, e.form = -0.25,
                  mechanism = c("E", "EC", "CE"),
                  ko = 1, kcf = 0, kcr = 0,
                  n = 1, alpha = 0.50, d = 1e-5, area = 0.01, 
                  temp = 298.15, scan.rate = 1.0, conc.bulk = 1e-3,
                  t.units = 1000, x.units = 100, sd.noise = 0, 
                  stir.rate = 0) {
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }
  
  mechanism = match.arg(mechanism)

  if (stir.rate > 10 | stir.rate < 0) {
    stop("Stir rate is a value between 0 and 10.")
  }
  
  if (mechanism == "E") {
    if (kcf != 0 | kcr != 0) {
      stop("For the E mechanism, kcf and kcr must have values of 0.")
    }
  }
  
  if (mechanism == "CE") {
    if (kcf <= 0) {
      stop("For the CE mechanism, kcf must have a value greater than 0.")
    }
  }
  
#' translate stir.rate to a factor that limits the thickness of the diffusion layer when stirring is turned on; this is accomplished in the loops for calculating the diffusion grid by dividing x.units by the scaling.factor defined here 
  
  scaling.factor = stir.rate * 5 + 1
  # if (stir.rate == "off") {
  #   scaling.factor = 1
  # } else if (stir.rate == "slow") {
  #   scaling.factor = 5
  # } else if (stir.rate == "medium") {
  #   scaling.factor = 15
  # } else if (stir.rate == "fast") {
  #   scaling.factor = 30
  # } else {
  #   stop("Options for stir rate are off, slow, medium, or fast.")
  # }
  
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
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = conc.bulk - cchem.bulk
      cred.bulk = 0
    } else {
      cox.bulk = conc.bulk
      cred.bulk = 0
      cchem.bulk = 0
    }
  } else {
    direction = +1
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = 0
      cred.bulk = conc.bulk - cchem.bulk
    } else {
      cox.bulk = 0
      cred.bulk = conc.bulk
      cchem.bulk = 0
    }
  }
  
  if (mechanism != "E") {
    min_tunits = 4 * t.tot * kcf
    if (t.units < min_tunits) {
      stop(paste("Minimum time units is", min_tunits, "if kcf =", kcf, "and with a total scan time of", t.tot, "s."))
    }
  }
  
  potential = seq(e.start, e.end, direction * scan.rate * delta.t) 
  
  kf = ko * exp(-alpha * n * f * (potential - e.form)/(r*temp))
  kb = ko * exp((1 - alpha) * n * f * (potential - e.form)/(r*temp))
  
  dif.ox = matrix(cox.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.red = matrix(cred.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  dif.chem = matrix(cchem.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
  if (mechanism == "CE") {
    for (i in 2:(t.units + 1)){
      for (j in 2:(x.units/scaling.factor)) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) + kcf * delta.t * dif.chem[i-1, j] - kcr * delta.t * dif.ox[i-1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.ox[i-1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) 
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) + kcf * delta.t*dif.chem[i-1, j] - kcr * delta.t * dif.red[i-1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.red[i-1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  } else {
    for (i in 2:(t.units + 1)){
      for (j in 2:(x.units/scaling.factor)) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kcf * delta.t * dif.red[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.red[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kcf * delta.t * dif.ox[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.ox[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  }
  
  noise = rnorm(t.units + 1, mean = 0,
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
  output = list("expt" = "LSV",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox * 10^6, 
                "reddata" = dif.red * 10^6, 
                "chemdata" = dif.chem * 10^6, 
                "formalE" = e.form,
                "initialE" = e.start,
                "endE" = e.end,
                "electrons" = n,
                "ko" = ko,
                "kcf" = kcf,
                "kcr" = kcr,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "scanrate" = scan.rate,
                "conc.bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction,
                "stir_rate" = stir.rate,
                "k_f" = kf,
                "k_b" = kb,
                "jox" = jox,
                "jred" = jred
  )
  invisible(output)
}

#' caSim: function to simulate a chronoamperometry experiment; parameters are similar to those in the functions cvSim with the addition of the parameter pulses, which takes a value of 1 for a single pulse experiment and a value of 2 for a double pulse experiment; the total elapsed time to the first pulse is given by t.1, the total elapsed time to the second pulse, if used, is t.2, and the total elapsed time to the experiment's end is t.end; see comments for the function cvSim for additional details not noted here

caSim = function(e.start = 0.0, e.pulse = -0.5, e.form = -0.25,
                 mechanism = c("E", "EC", "CE"),
                 ko = 1, kcf = 0, kcr = 0,
                 pulses = c("single", "double"), 
                 t.1 = 10, t.2 = 0, t.end = 30,
                 n = 1, alpha = 0.50, d = 1e-5, area = 0.01, 
                 temp = 298.15, scan.rate = 1.0, conc.bulk = 1e-3,
                 t.units = 1000, x.units = 100, sd.noise = 0) {
  
  if (x.units >= sqrt(18 * t.units)) {
    stop("x.units must be less than sqrt(18 * t.units)")
  }
  
  mechanism = match.arg(mechanism)
  pulses = match.arg(pulses)
  
  if (mechanism == "E") {
    if (kcf != 0 | kcr != 0) {
      stop("For the E mechanism, kcf and kcr must have values of 0.")
    }
  }
  
  if (mechanism == "CE") {
    if (kcf <= 0) {
      stop("For the CE mechanism, kcf must have a value greater than 0.")
    }
  }
  
  if (pulses == "single") {
    num.pulses = 1
    if (t.2 != 0) {
      stop("For a single pulse experiment, t.2 must be 0.")
    }
  } else {
    num.pulses = 2
    if (t.2 <= t.1) {
      stop("For a double pulse experiment, t.2 must be greater than t.1")
    }
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
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = conc.bulk - cchem.bulk
      cred.bulk = 0
    } else {
      cox.bulk = conc.bulk
      cred.bulk = 0
      cchem.bulk = 0
    }
  } else {
    direction = +1
    if (mechanism == "CE") {
      cchem.bulk = conc.bulk/(1 + kcf/kcr)
      cox.bulk = 0
      cred.bulk = conc.bulk - cchem.bulk
    } else {
      cox.bulk = 0
      cred.bulk = conc.bulk
      cchem.bulk = 0
    }
  }
  
  if (mechanism != "E") {
    min_tunits = 4 * num.pulses * (kcf + kcr)
    # min_tunits = 4 * t.tot * kcf
    if (t.units < min_tunits) {
      stop(paste("Minimum time units is", min_tunits, "for these conditions."))
    }
  }
  
  if (pulses == "single") { 
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
  dif.chem = matrix(cchem.bulk/1000, nrow = t.units + 1, ncol = x.units + 1)
  
  jox = rep(0, t.units + 1)
  jred = rep(0, t.units + 1)
  current.total = rep(0, t.units + 1)
  
  if (mechanism == "CE") {
    for (i in 2:(t.units + 1)){
      for (j in 2:x.units) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) + kcf * delta.t * dif.chem[i-1, j] - kcr * delta.t * dif.ox[i-1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.ox[i-1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) 
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) + kcf * delta.t*dif.chem[i-1, j] - kcr * delta.t * dif.red[i-1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) - kcf * delta.t * dif.chem[i-1, j] + kcr * delta.t * dif.red[i-1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  } else {
    for (i in 2:(t.units + 1)){
      for (j in 2:x.units) {
        if (direction == -1) {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1])
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) - kcf * delta.t * dif.red[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.red[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        } else {
          dif.ox[i, j] = dif.ox[i-1, j] + lambda * (dif.ox[i-1, j-1] - 2 * dif.ox[i-1, j] + dif.ox[i-1, j+1]) - kcf * delta.t * dif.ox[i-1, j] + kcr * delta.t * dif.chem[i - 1, j]
          dif.red[i, j] = dif.red[i-1, j] + lambda * (dif.red[i-1, j-1] - 2 * dif.red[i-1, j] + dif.red[i-1, j+1]) 
          dif.chem[i, j] = dif.chem[i-1, j] + lambda * (dif.chem[i-1, j-1] - 2 * dif.chem[i-1, j] + dif.chem[i-1, j+1]) + kcf * delta.t * dif.ox[i-1, j] - kcr * delta.t * dif.chem[i - 1, j]
        }
      }
      jox[i] = -(kf[i] * dif.ox[i,2] - kb[i] * dif.red[i,2])/(1 + (kf[i] * delta.x)/d + (kb[i] * delta.x)/d)
      jred[i] = -jox[i]
      dif.ox[i, 1] = dif.ox[i, 2] + jox[i] * delta.x/d
      dif.red[i, 1] = dif.red[i, 2] + jred[i] * delta.x/d
      if (direction == -1) {
        dif.chem[i,1] = dif.chem[i, 2]
      } else {
        dif.chem[i,1] = dif.chem[i, 2]
      }
      current.total[i] = -n * f * area * jox[i]
    }
  }
  
  noise = rnorm(t.units + 1, mean = 0,
                sd = sd.noise * max(abs(current.total))/100)
  current.total = current.total + noise
  
  output = list("expt" = "CA",
                "mechanism" = mechanism,
                "file_type" = "full",
                "current" = current.total*10^6, 
                "potential" = potential,
                "time" = time, 
                "distance" = distance, 
                "oxdata" = dif.ox * 10^6, 
                "reddata" = dif.red * 10^6, 
                "chemdata" = dif.chem * 10^6,
                "formalE" = e.form,
                "initialE" = e.start,
                "pulseE" = e.pulse,
                "electrons" = n,
                "ko" = ko,
                "kcf" = kcf,
                "kcr" = kcr,
                "alpha" = alpha,
                "diffcoef" = d,
                "area" = area,
                "temperature" = temp,
                "conc.bulk" = conc.bulk,
                "tunits" = t.units,
                "xunits" = x.units,
                "sdnoise" = sd.noise,
                "direction" = direction,
                "pulses" = pulses,
                "time_pulse1" = t.1,
                "time_pulse2" = t.2,
                "time_end" = t.end,
                "k_f" = kf,
                "k_b" = kb,
                "jox" = jox,
                "jred" = jred
  )
  invisible(output)
}

#' ccSim: function that converts the result of a chronoamperometry simulation into its corresponding chonocoulometry simulation by integrating current over time using trapazoidal rule integration (adapted here from the function trapz in the caTools package)

ccSim = function(filename){
  
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
                "file_type" = filename$file_type,
                "charge" = charge, 
                "potential" = filename$potential,
                "time" = filename$time, 
                "distance" = filename$distance, 
                "oxdata" = filename$oxdata, 
                "reddata" = filename$reddata,
                "chemdata" = filename$chemdata,
                "formalE" = filename$formalE,
                "initialE" = filename$initialE,
                "pulseE" = filename$pulseE,
                "electrons" = filename$electrons,
                "ko" = filename$ko,
                "kcf" = filename$kcf,
                "kcr" = filename$kcr,
                "alpha" = filename$alpha,
                "diffcoef" = filename$diffcoef,
                "area" = filename$area,
                "temperature" = filename$temperature,
                "conc.bulk" = filename$conc.bulk,
                "tunits" = filename$tunits,
                "xunits" = filename$xunits,
                "sdnoise" = filename$sdnoise,
                "direction" = filename$direction,
                "pulses" = filename$pulses,
                "time_pulse1" = filename$time_pulse1,
                "time_pulse2" = filename$time_pulse2,
                "time_end" = filename$time_end,
                "k_f" = filename$kf,
                "k_b" = filename$kb,
                "jox" = filename$jox,
                "jred" = filename$jred
  )
  
  invisible(output)
  
}

#' sampleVoltgram: used to create a reduced data file of potentials and currents for a cyclic voltammogram or a linear sweep voltammogram; data.reduction gives the percentage of data to keep, which then is evenly spaced using a calculated increment 

sampleVoltgram = function(filename, data.reduction = 1) {
  
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
                  "ko" = filename$ko,
                  "kcf" = filename$kcf,
                  "kcr" = filename$kcr,
                  "alpha" = filename$alpha,
                  "diffcoef" = filename$diffcoef,
                  "area" = filename$area,
                  "temperature" = filename$temperature,
                  "scanrate" = filename$scanrate,
                  "conc.bulk" = filename$conc.bulk,
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
                  "ko" = filename$ko,
                  "kcf" = filename$kcf,
                  "kcr" = filename$kcr,
                  "alpha" = filename$alpha,
                  "diffcoef" = filename$diffcoef,
                  "area" = filename$area,
                  "temperature" = filename$temperature,
                  "scanrate" = filename$scanrate,
                  "conc.bulk" = filename$conc.bulk,
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

#' sampleAmpgram: used to create a reduced data file of times and currents for a chronoamperogram; data.reduction gives the percentage of data to keep, which then is evenly spaced using a calculated increment 

sampleAmpgram = function(filename, data.reduction = 1) {
  
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
                "ko" = filename$ko,
                "kcf" = filename$kcf,
                "kcr" = filename$kcr,
                "alpha" = filename$alpha,
                "diffcoef" = filename$diffcoef,
                "area" = filename$area,
                "temperature" = filename$temperature,
                "conc.bulk" = filename$conc.bulk,
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

#' sampleCoulgram: used to create a reduced data file of times and charges for a chronocoulogram; data.reduction gives the percentage of data to keep, which then is evenly spaced using a calculated increment 

sampleCoulgram = function(filename, data.reduction = 1) {
  
  time = filename$time
  charge = filename$charge
  len = length(filename$time)
  delta.data = len * data.reduction/100
  time.new = time[seq(1, len, len/delta.data)]
  charge.new = charge[seq(1, len, len/delta.data)]
  output = list("expt" = filename$expt,
                "mechanism" = filename$mechanism,
                "file_type" = "reduced",
                "charge" = charge.new, 
                "time" = time.new,
                "formalE" = filename$formalE,
                "initialE" = filename$initialE,
                "pulseE" = filename$pulseE,
                "electrons" = filename$electrons,
                "ko" = filename$ko,
                "kcf" = filename$kcf,
                "kcr" = filename$kcr,
                "alpha" = filename$alpha,
                "diffcoef" = filename$diffcoef,
                "area" = filename$area,
                "temperature" = filename$temperature,
                "conc.bulk" = filename$conc.bulk,
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
