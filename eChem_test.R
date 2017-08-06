#' Example 1: Gosser's test conditions for E mechanism; see box on page 36 and Figure 2-7 and Figure 2-8
example1 = cvSim(e.start = 0, e.switch = -0.5, e.form = -0.25, mechanism = "E", scan.rate = 1, area = 0.01, temp = 298, conc.bulk = 1e-3, n = 1, d = 1e-5, alpha = 0.5, ko = 1, kcf = 0, kcr = 0)
old.par = par(mfrow = c(2, 2))
plotPotential(example1)
plotVoltgram(list(example1))
plotDiffusion(example1, t = 0.25)
plotDiffusion(example1, t = 0.75)
par(old.par)
annotateCV(example1)
diffusionGrid(example1, species = "All", scale.factor = 5)

#' Example 2: Effect of the heterogeneous electron-transfer rate constant (ko) using conditions similar to those in Figure 2-11 in Gosser; note, Gosser does not give value of ko for reversible conditions, which is taken here as 1 cm/s
example2a = cvSim(ko = 1, e.switch = -0.8, e.form = -0.4)
example2b = cvSim(ko = 0.01, e.switch = -0.8, e.form = -0.4)
example2c = cvSim(ko = 0.001, e.switch = -0.8, e.form = -0.4)
example2d = cvSim(ko = 0.0001, e.switch = -0.8, e.form = -0.4)
plotVoltgram(list(example2a, example2b, example2c, example2d), legend_text = c("ko = 1 cm/s", "ko = 0.01 cm/s", "ko = 0.001 cm/s", "ko = 0.0001 cm/s"))

#' Example 3: Effect of transfer coefficient (alpha) using conditions similar to those in Figure 2-12 in Gosser
example3a = cvSim(ko = 0.0001, alpha = 0.75, e.start = 0.4, e.switch = -1.2, e.form = -0.4)
example3b = cvSim(ko = 0.0001, alpha = 0.50, e.start = 0.4, e.switch = -1.2, e.form = -0.4)
example3c = cvSim(ko = 0.0001, alpha = 0.25, e.start = 0.4, e.switch = -1.2, e.form = -0.4)
plotVoltgram(list(example3a, example3b, example3c), legend_text = c("alpha = 0.75", "alpha = 0.50", "alpha = 0.25"))

#' Example 4: Gosser's test conditions for EC mechanism; see box on page 48 and Figure 2-14; for the last example it is necessary to increase the number of time units
example4a = cvSim(mechanism = "EC", ko = 1, kcf = 0, kcr = 0)
example4b = cvSim(mechanism = "EC", ko = 1, kcf = 2.5, kcr = 0)
example4c = cvSim(mechanism = "EC", ko = 1, kcf = 100, kcr = 0)
example4d = cvSim(mechanism = "EC", ko = 1, kcf = 1000, kcr = 0, t.units = 4000)
plotVoltgram(list(example4a, example4b, example4c, example4d), legend_text = c("kcf = 0", "kcf = 2.5", "kcf = 100", "kcf = 1000"))

#' Example 5: Effect of kcf on results of chronoamperometry experiments; see Gosser's Figure 2-26; note that the times here are offset by 1 s relative to figure
example5a = caSim(mechanism = "EC", pulses = "double", t.1 = 1, t.2 = 3, t.end = 5, kcf = 50)
example5b = caSim(mechanism = "EC", pulses = "double",t.1 = 1, t.2 = 3, t.end = 5, kcf = 100)
example5c = caSim(mechanism = "EC", pulses = "double",t.1 = 1, t.2 = 3, t.end = 5, kcf = 500)
plotAmpgram(list(example5a, example5b, example5c))

#' Example 6: Duplication of Figure 2.12 in Brownson, D. A. C.; Banks, C. E. The Handbook of Graphene Electrochemistry; note that the CV here maintains an x-axis that runs from more +E to more -E and that the individual diffusion profiles are at different times than in Figure 2.12
example6 = cvSim(e.start = -0.4, e.switch = +0.4, e.form = 0, area = 1)
plotGrid(example6)

#' Example 7: Verification of square root dependence of current on scan rate for CV with E mechanism
example7a = cvSim(scan.rate = 0.01)
example7b = cvSim(scan.rate = 0.05)
example7c = cvSim(scan.rate = 0.1)
example7d = cvSim(scan.rate = 0.5)
example7e = cvSim(scan.rate = 1)
example7f = cvSim(scan.rate = 5)
example7g = cvSim(scan.rate = 10)
plotVoltgram(list(example7a, example7c, example7e, example7g), legend_text = c("v = 0.01", "v = 0.1", "v = 1", "v = 10"))
scanrates = c(0.01, 0.05, 0.1, 0.5, 1, 5, 10)
peakcurrents = c(max(example7a$current), max(example7b$current),max(example7c$current), max(example7d$current), max(example7e$current), max(example7f$current), max(example7g$current))
plot(sqrt(scanrates), peakcurrents, pch = 19, col = "blue", xlab = "(scan rate)^0.5", ylab = "peak current (µA)")
iv.lm = lm(peakcurrents ~ sqrt(scanrates))
abline(iv.lm, col = "blue", lty = 1)

#' Example 8: Duplicating results shown in left half of Figure 24 of Rick Kelly's ASDL module
example8a = cvSim(mechanism = "CE", e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 0.1, ko = 100, kcf = 0.1, kcr = 1, area = 0.10)
example8b = cvSim(mechanism = "CE",e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 0.1, ko = 100, kcf = 0.1, kcr = 0.1, area = 0.10)
example8c = cvSim(mechanism = "CE",e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 0.1, ko = 100, kcf = 0.1, kcr =0.01, area = 0.10)
plotVoltgram(list(example8a, example8b, example8c), legend_text = c("K = 0.1", "K = 1", "K = 10"))

#' Example 9: Duplicating results shown in right half of Figure 24 of Rick Kelly's ASDL module; note that Figure 24 probably mislabels K as 10 instead of 0.10
example9a = cvSim(mechanism = "CE",e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 0.01, ko = 100, kcf = 0.1, kcr = 1, area = 0.10)
example9b = cvSim(mechanism = "CE",e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 0.1, ko = 100, kcf = 0.1, kcr = 1, area = 0.10)
example9c = cvSim(mechanism = "CE",e.start = 0, e.switch = -0.7, e.form = -0.3, scan.rate = 1, ko = 100, kcf = 0.1, kcr = 1, area = 0.10)
plotVoltgram(list(example9a, example9b, example9c), legend_text = c("v = 0.01 V/s", "v = 0.1 V/s", "v = 1 V/s"))
