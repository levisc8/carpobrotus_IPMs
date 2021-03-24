# Diagnostics for shrinkage models

# precip plots
hs_diagnostic_plots("prec", "repro", "lin") 
hs_diagnostic_plots("prec", "repro", "quad") 
hs_diagnostic_plots("prec", "flower_n", "lin") 
hs_diagnostic_plots("prec", "flower_n", "quad") 


# tmax
hs_diagnostic_plots("tmax", "repro", "lin") 
hs_diagnostic_plots("tmax", "repro", "quad") 
hs_diagnostic_plots("tmax", "flower_n", "lin") 
hs_diagnostic_plots("tmax", "flower_n", "quad") 


# tmin
hs_diagnostic_plots("tmin", "repro", "lin") 
hs_diagnostic_plots("tmin", "repro", "quad") 
hs_diagnostic_plots("tmin", "flower_n", "quad") 
hs_diagnostic_plots("tmin", "flower_n", "lin") 
