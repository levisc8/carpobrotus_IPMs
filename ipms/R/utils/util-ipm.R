# rename fixed effects from brms to ipmr names

name_fixed_pars <- function(model, vr, sw = 3) {
  
  mod_nm <- deparse(substitute(model))
  
  base_nms <- paste(vr, 
                    c(
                      "int", "z",
                      "temp_dry", "temp_wet", 
                      "prec_dry", "prec_wet",
                      paste(paste0("sw", sw), 
                            c("dry", "wet"), sep = "_"),
                      "native"
                    ), sep = "_"
  )
  
  interaction_nms <- rownames(fixef(model)) %>%
    .[grepl("\\:", .)] %>%
    gsub(pattern = "\\:", replacement = "_x_", x = .) %>%
    gsub(pattern = "_t$|_t_1$", replacement = "", x = .)
  interaction_nms <- paste(vr, interaction_nms, sep = "_")
  
  out <- c(base_nms, interaction_nms)
  
  if(mod_nm == "grow_mod") {
    
    base_nms <- c(base_nms[1],
                  paste0(vr, "_sigma_int"), 
                  base_nms[2:length(base_nms)])
    
    sig_nms <- rownames(fixef(model)) %>%
      .[grepl("sigma_", .) & !grepl("Intercept", .)] %>%
      gsub(pattern = "_t$$", replacement = "", x = .)
    sig_nms <- paste0("grow_", sig_nms)
    
    out <- c(base_nms, interaction_nms, sig_nms)
  }
  
  out <- setNames(fixef(model)[ , 1], out)
  
  return(as.list(out))
  
}

# rename random effects from brms to ipmr names

name_ran_pars <- function(model, vr) {
  
  mod_nm <- deparse(substitute(model))
  
  pops <- rownames(ranef(model)[[1]][ , , 1])
  
  base_nm <- paste(vr, "intercept", pops, sep = "_")
  
  out <- setNames(ranef(model)[[1]][ , 1, 1], base_nm)
  
  if(mod_nm == "grow_mod") {
    sig_nms <- rownames(ranef(model)[[1]][ , , 2])
    sig_nms <- paste0("grow_sigma_intercept_", sig_nms)
    sig_out <- setNames(ranef(model)[[1]][ , 1, 2], sig_nms)
    
    out <- c(out, sig_out)
  }
  
  return(as.list(out))
  
}

# Rename posterior draws data.frames to ipmr

rename_draws <- function(draws, vr) {
  
  nm <- deparse(substitute(draws))
  
  # Fix random effects names
  names(draws)[grepl("\\,Intercept\\]", names(draws))] <- 
    gsub("\\,Intercept\\]", 
         "", 
         names(draws)[grepl("\\,Intercept\\]", names(draws))])
  
  names(draws)[grepl("r_site\\[", names(draws))] <-
    gsub("r_site\\[",
         paste0(vr, "_intercept_"), 
         names(draws)[grepl("r_site\\[", names(draws))])
  
  # Site random effect for growth variance model
  names(draws)[grepl("r_site__sigma\\[", names(draws))] <-
    gsub("r_site__sigma\\[",
         paste0(vr, "_sigma_intercept_"), 
         names(draws)[grepl("r_site", names(draws))])
  
  # Fix the fixed effects now.
  names(draws)[!grepl(vr, names(draws))] <- 
    paste(vr,
          names(draws)[!grepl(vr, names(draws))], sep = "_")
  
  # Remove leading 'b's (now no longer leading)
  
  names(draws) <- gsub("_b_", "_", names(draws))
  
  # Interaction terms
  names(draws) <- gsub("\\:", "_x_", names(draws))
  
  # climate terms
  names(draws) <- gsub("_t$|_t_1$", "", names(draws))
  
  # Fix intercepts
  
  names(draws) <- gsub("Intercept", "int", names(draws))
  # Fix the name of the size coefficient
  z_nm <- paste0(vr, "_log_size")
  
  names(draws)[names(draws) == z_nm] <- paste0(vr, "_z")
  
  
  return(draws)
}

# Simple switch function that takes a symbol site name and returns native status

is_native <- function(x) {
  x <- s_char(x)
  if(x %in% c("Struisbaai", "Rooisand",
              "Vogelgat", "Springfontein",
              "Melkboss", "St_Francis")) {
    return(1)
  } else {
    return(0)
  }
}

# Symbol -> character()

s_char <- function(x) {
  enexpr(x) %>% expr_text()
}

# Returns predicted values w/o CIs from brms

predict_carp_ipm <- function(model, data) {
  x <- predict(model, newdata = data, cores = 4L)
  unname(x[ , 1, drop = TRUE])
}

# Either summarizes the draws to mean posterior value, or selects a specific
# draw (e.g. for uncertainty simulation)
prep_draws <- function(par, draws, draw_id, f = mean) {
  
  if(!is.null(draw_id)){
    
    draws[grepl(par, names(draws))] %>%
      .[draw_id, ] %>%
      unlist()
    
  } else {
    
    draws[grepl(par, names(draws))] %>%
      apply(2, f)
    
  }
}

