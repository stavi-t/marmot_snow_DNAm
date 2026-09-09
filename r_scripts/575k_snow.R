library(lme4breeding)
library(parallel)
library(nadiv)
library(ybamaRmot)
library(abind)
library(tidyverse)

met_dat <- read_csv("metadata_springT_189.csv") %>%
  mutate(
    uid = as.factor(uid),
    col_area = as.factor(col_area),
    year = as.factor(year),
    uidp = uid,
    batch = as.factor(batch),
    site = file_name,
    obs = seq_len(nrow(.))
    #   file_name = gsub("-", ".", file_name),
  )

################################################################################
## can be skipped (slow) and go straight to read RDS below
## loading  sites data
sites <- read.delim("raw_data/575k_0.80_165_count.txt", sep = " ", row.names = 1)
# sites <- sites %>%
#  mutate(
#    Site = gsub("-", ".", Site),
#    Site = gsub(" ", "_", Site)
#  )

tot_sites <- nrow(sites)


## estimating site variance and mu
meth <- grep("meth", colnames(sites))
tot <- grep("total", colnames(sites))
site_var <- data.frame(site = c(sites$Site))
site_var$nb <- unlist(
  mclapply(
    1:tot_sites,
    function(x) {
      sum(!is.na(unique(
        unlist(sites[x, meth]) /
          unlist(sites[x, tot])
      )))
    },
    mc.cores = 8
  )
)
write_csv(site_var, "575k_output/site_variance.csv")

site_var$var <- unlist(
  mclapply(
    1:tot_sites,
    function(x) {
      var(
        unlist(sites[x, meth]) /
          unlist(sites[x, tot]),
        na.rm = TRUE
      )
    },
    mc.cores = 8
  )
)
write_csv(site_var, "575k_output/site_variance.csv")

site_var$n_na <- unlist(
  mclapply(
    1:tot_sites,
    function(x) {
      sum(is.na(unlist(sites[x, meth])))
    },
    mc.cores = 8
  )
)
write_csv(site_var, "575k_output/site_variance.csv")

site_var$mu <- unlist(
  mclapply(
    1:tot_sites,
    function(x) {
      mean(
        unlist(sites[x, meth]) /
          unlist(sites[x, tot]),
        na.rm = TRUE
      )
    },
    mc.cores = 8
  )
)
write_csv(site_var, "575k_output/site_variance.csv")

sites <- inner_join(sites, site_var, by = join_by(Site == site))

saveRDS(sites, "575k_output/sites.rds")

################################################################################

sites <- readRDS("575k_output/sites.rds")

site_sel <- filter(sites, nb > 30 & var > 0.002)
site_names <- c(site_sel$Site)
rm(sites)

ped <- prepPed(fn_ped())
ped <- nadiv::prunePed(ped, phenotyped = met_dat$uid)

A <- makeA(ped[, 1:3])

run_model <- function(x, param) {
  tots <- site_sel[x, ] %>%
    select(contains("total")) %>%
    pivot_longer(
      cols = contains("total"),
      names_to = "site",
      values_to = "tot"
    ) %>%
    mutate(site = gsub("_total_cytosine", "", site))
  cs <- site_sel[x, ] %>%
    select(contains("meth")) %>%
    pivot_longer(
      cols = contains("meth"),
      names_to = "site",
      values_to = "cs"
    ) %>%
    mutate(site = gsub("_meth_cytosine", "", site))
  dat <- full_join(tots, cs, by = "site") %>%
    inner_join(met_dat, by = join_by(site)) %>%
    mutate(
      ts = tot - cs
    )
  assign("dat_x", dat, envir = .GlobalEnv)
  if (param == "snw") {
    assign(
      "form",
      as.formula("cbind(cs, ts) ~ age + snow_melt + batch + (1 | uid) + (1 | uidp) + (1 | obs)"),
      envir = .GlobalEnv
    )
  } else if (param == "sprT") {
    #bug in lme4breed can't read from function environment
        assign(
          "form",
          as.formula("cbind(cs, ts) ~ age + spr_avg_temp + batch + (1 | uid) + (1 | uidp) + (1 | obs)"),
          envir = .GlobalEnv
        )
   } else {
    stop("param must be snw or sprT")
  }

  m <- tryCatch(
    {
      m0 <- lmebreed(
        form,
        relmat = list(uid = A),
        data = dat_x,
        family = "binomial",
        verbose = FALSE,
        dateWarning = FALSE
      )
      summary(m0)
    },
    error = function(war) {
      return(NA)
    }
  )
  m
}



fn_site_mc <- function(site_num, n_cores, params) {
  site_num <- site_num[site_num <= nrow(site_sel)]
  out <- mclapply(
    site_num,
    function(x) run_model(x, params),
    mc.cores = n_cores
  )
  names(out) <- site_sel$Site[site_num]
  out
}
by_site <- fn_site_mc(1:16, 8, params = "snw")


fn_site <- function(site_num) {
  site_num <- site_num[site_num <= nrow(site_sel)]
  out <- lapply(
    site_num,
    run_model
  )
  names(out) <- site_sel$Site[site_num]
  out
}

for (i in 1:6) {
  print(i)
  by_site <- fn_site_mc(((i - 1) * 100000 + 1):(i * 100000), 6, "snw")
  saveRDS(by_site, paste0("575k_output/lmebreed_site_snow_m", i, ".rds"))
  rm(by_site)
}

for (i in 1:6) {
  print(i)
  by_site <- fn_site_mc(((i - 1) * 100000 + 1):(i * 100000), 6, "sprT")
  saveRDS(by_site, paste0("575k_output/lmebreed_site_sprT_m", i, ".rds"))
  rm(by_site)
}

###############################################
# extracting estimates


p_age <- lapply(by_site, function(x) {
  x$coefficients[2, 4]
})
hist(unlist(p_age))

n <- dimnames(by_site[[1]]$coefficients)

for (i in 1:6) {
  by_site <- readRDS(paste0("575k_output/lmebreed_site_sprT_m", i, ".rds"))
  out_array <- array(
    unlist(lapply(
      by_site, function(x) {
        x$coefficients
      }
    )),
    c(5, 4, length(by_site)),
    dimnames = list(n[[1]], n[[2]], names(by_site))
  )
  if (i == 1) {
    coefs_all_sites <- out_array
  } else {
    coefs_all_sites <- abind(coefs_all_sites, out_array, along = 3)
  }
}
rm(by_site, out_array)

output <- aperm(coefs_all_sites, c(3, 2, 1))
saveRDS(output, "575k_output/fixed_coefs_575k_sprT.rds")




################################################################################
# playing with the results
output <- readRDS("575k_output/fixed_coefs_575k_snow.rds")

## in array notation [x,y,z] x is a mthylation site, y is type of estimate (estimate, se, t p) and z is the parameter, 3 being snowmelt and 2 age


## for snowmelt
hist(output[, 4, 3])
sum(output[, 4, 3] < 0.001)
