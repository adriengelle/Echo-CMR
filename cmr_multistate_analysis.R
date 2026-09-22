###############################################################################
## Echo parakeet - Multistate CMR analysis (RMark / Program MARK)
##
## Fits an Arnason-Schwarz multistate Cormack-Jolly-Seber (CJS) model to estimate
## age- and time-dependent survival, taking into account heterogeneity in 
## detection and transition probabilities, then produces a juvenile 
## (<2 years old) / adult survival figure by population and a model summary.
##
## REQUIREMENTS: R, Program MARK installed, and the packages loaded below.
## See 'README.md' for full setup and step-by-step instructions.
##
## HOW TO RUN: Source the whole file (Ctrl+Shift+S) or click "Source" after 
## adding encounter history file to the '/data' folder
## (script and data folder should be on the root of working directory)
###############################################################################

library(RMark) #interface to Program MARK  (requires MARK installed)
library(dplyr)
library(ggplot2)

## ============================ CONFIG ========================================
## Edit these if necessary, then source the script.
output_dir <- "output"
input_dir <- "data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) #create directories if non-existing
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
file <- list.files(path = "data/", pattern = "\\.csv$", full.names = TRUE)[1]
data_file <- read.csv(file) #read unique file in 'data'

## State codes used in the year columns of the CSV -> numeric strata for MARK
## PB = pre-breeder (1), B = breeder (2), PO = post-breeder (3), 0 = not seen
state_codes <- c("B" = "2", "PB" = "1", "PO" = "3", "0" = "0")
## ===========================================================================

## ---- 1. Read and prepare the encounter history ----------------------------
EH <- read.csv(file, header = TRUE, check.names = FALSE)
EH[] <- lapply(EH, function(x) if (is.character(x)) trimws(x) else x) #removes leading/trailing whitespace (state attribution)

if (!"chrt" %in% names(EH))
  stop("The CSV has no 'chrt' (cohort) column, which is needed to set the years.")
EH$chrt <- as.numeric(as.character(EH$chrt))
max_chrt         <- max(EH$chrt, na.rm = TRUE)
start_year       <- min(EH$chrt, na.rm = TRUE)  #first survey year = earliest ringed cohort
end_year         <- max_chrt + 1                #last survey year  = year after last cohort
last_cohort_year <- max_chrt - 1                #analyse cohorts up to one year before the last

if ("subpop" %in% names(EH)) {
  EH$subpop <- as.factor(EH$subpop)
  EH <- subset(EH, subpop %in% c("BO", "GG"))
}

ch_columns <- as.character(start_year:end_year)
if (!all(ch_columns %in% names(EH)))
  stop("CSV is missing year columns: ",
       paste(setdiff(ch_columns, names(EH)), collapse = ", "),
       "\nCheck start_year / end_year against your file header.")

ch <- apply(EH[, ch_columns, drop = FALSE], 1, function(x) {
  codes <- state_codes[as.character(x)]
  codes[is.na(codes)] <- "0"          #anything unrecognised -> not seen
  paste(codes, collapse = "")
})

rmark_data <- data.frame(
  ch = ch,
  EH[, !names(EH) %in% c(ch_columns, "id"), drop = FALSE],
  stringsAsFactors = FALSE
)

all_zero <- which(rmark_data$ch == paste(rep("0", length(ch_columns)), collapse = ""))
if (length(all_zero) > 0) rmark_data <- rmark_data[-all_zero, ]
if (nrow(rmark_data) == 0) stop("No usable encounter histories after filtering.")

## ---- 2. Process data + build design data (Multistrata) --------

ms.proc <- process.data(rmark_data, model = "Multistrata",
                        strata.labels = c("1", "2", "3"), group = 'subpop')#by pop
ms.ddl  <- make.design.data(ms.proc)

## ---- 3. Age classes and biologically impossible cells (fixed to 0) ---------
{
# Survival (S)
ms.ddl$S$ageclass3 <- cut(ms.ddl$S$Age, c(0, 2, 4, Inf), right = FALSE)
levels(ms.ddl$S$ageclass3) <- c("Juveniles(0-1)", "Adults(2-3)", "Adults(4+)")
ms.ddl$S$juv         <- ifelse(ms.ddl$S$ageclass3 == "Juveniles(0-1)", 1, 0)
ms.ddl$S$adults      <- ifelse(ms.ddl$S$Age >= 2, 1, 0)

ms.ddl$S$fix <- NA
ms.ddl$S$fix[ms.ddl$S$Age == 1 & ms.ddl$S$stratum != "1"] <- 0
ms.ddl$S$fix[ms.ddl$S$Age == 2 & ms.ddl$S$stratum == "3"] <- 0

# Detection (p)
ms.ddl$p$ageclass2 <- cut(ms.ddl$p$Age, c(-Inf, 2, Inf), right = FALSE)
levels(ms.ddl$p$ageclass2) <- c("Young", "Age2plus")
ms.ddl$p$age1     <- ifelse(ms.ddl$p$ageclass2 == "Young", 1, 0)
ms.ddl$p$age2plus <- ifelse(ms.ddl$p$ageclass2 == "Age2plus", 1, 0)

ms.ddl$p$fix <- NA
ms.ddl$p$fix[ms.ddl$p$stratum == "2" & ms.ddl$p$Age == 1] <- 0
ms.ddl$p$fix[ms.ddl$p$stratum == "3" & ms.ddl$p$Age < 3] <- 0

# Transitions (Psi)
ms.ddl$Psi$age1     <- ifelse(ms.ddl$Psi$Age == 1, 1, 0)
ms.ddl$Psi$age2     <- ifelse(ms.ddl$Psi$Age == 2, 1, 0)
ms.ddl$Psi$age3plus <- ifelse(ms.ddl$Psi$Age >= 3, 1, 0)
ms.ddl$Psi$growth_period <- cut(ms.ddl$Psi$Time, breaks = c(0, 10, 22, Inf),
                                labels = c("low", "steep_incline", "plateau"),
                                right = FALSE)

s  <- as.character(ms.ddl$Psi$stratum)
ts <- as.character(ms.ddl$Psi$tostratum)
ms.ddl$Psi$trans_1to2 <- as.integer(s == "1" & ts == "2")
ms.ddl$Psi$trans_2to3 <- as.integer(s == "2" & ts == "3")
ms.ddl$Psi$trans_3to2 <- as.integer(s == "3" & ts == "2")

#impossible transitions (1 to 3, 2 to 1, 3 to 1)
ms.ddl$Psi$fix <- NA
ms.ddl$Psi$fix[s == "1" & ts == "3"] <- 0
ms.ddl$Psi$fix[s == "2" & ts == "1"] <- 0
ms.ddl$Psi$fix[s == "3" & ts == "1"] <- 0
ms.ddl$Psi$fix[s == "1" & ts == "2" & ms.ddl$Psi$Age < 1] <- 0 #can't be breeder before age 2
ms.ddl$Psi$fix[s == "2" & ts == "3" & ms.ddl$Psi$Age < 2] <- 0
ms.ddl$Psi$fix[s == "3" & ts == "2" & ms.ddl$Psi$Age < 3] <- 0
}
## ---- 4. Fit the model (calls Program MARK) ----------------------------
ms_model <- mark(ms.proc, ms.ddl,
  model.parameters = list(
    S   = list(formula = ~ -1 + juv:time + adults:stratum:time + group),
    p   = list(formula = ~ age1:time + age2plus:stratum:time + group),
    Psi = list(formula = ~ -1 + trans_1to2:age1:growth_period +
                              trans_1to2:age2:growth_period +
                              trans_1to2:age3plus:growth_period +
                              trans_2to3:age2:growth_period +
                              trans_2to3:age3plus:growth_period +
                              trans_3to2:age3plus:growth_period)),
  output = FALSE, delete = TRUE)

## ---- 5. Model summary + real estimates -------------------------------------
sink(file.path(output_dir, "model_summary.txt"))
print(summary(ms_model))
cat("\nAICc:      ", ms_model$results$AICc, "\n")
cat("Deviance:  ", ms_model$results$deviance, "\n")
cat("Parameters:", ms_model$results$npar, "\n")
sink()

S.real   <- get.real(ms_model, "S",   se = TRUE)
p.real   <- get.real(ms_model, "p",   se = TRUE)
Psi.real <- get.real(ms_model, "Psi", se = TRUE)

#write in output folder
write.csv(S.real,   file.path(output_dir, "survival_estimates.csv"),   row.names = FALSE)
write.csv(p.real,   file.path(output_dir, "detection_estimates.csv"),  row.names = FALSE)
write.csv(Psi.real, file.path(output_dir, "transition_estimates.csv"), row.names = FALSE)

## ---- 6. Survival plot (juveniles 0-1 vs adults 2+) -------------------------
## Drop the last (confounded) estimable year per age class, 
# and transform 'occasion' to calendar year.
threshold <- 0.42 #juvenile survival threshold (PVA)

surv <- S.real %>%
  #filter(se > 0) %>%                                   # drop fixed / structural cells
  mutate(
    time    = as.numeric(as.character(time)),
    stratum = as.character(stratum),
    pop     = as.character(subpop),                    # check names(S.real) if this errors
    class   = case_when(
      Age <  2                   ~ "Juveniles (0-1 yrs)",
      Age >= 2 & stratum == "1"  ~ "Pre-breeders (2+ yrs)",
      Age >= 2 & stratum == "2"  ~ "Breeders (2+ yrs)",
      Age >= 2 & stratum == "3"  ~ "Post-breeders (3+ yrs)"
    )
  ) %>%
  filter(!is.na(class)) %>%
  distinct(class, pop, time, .keep_all = TRUE) %>%      # unique by class x pop x year
  group_by(class, pop) %>%
  filter(time < max(time)) %>%                          # drop last estimable year, per line
  ungroup() %>%
  mutate(
    Year  = start_year + time - 1,
    class = factor(class, levels = c("Juveniles (0-1 yrs)", "Pre-breeders (2+ yrs)",
                                     "Breeders (2+ yrs)", "Post-breeders (3+ yrs)")),
    Population = factor(ifelse(pop == "GG", "Gorges (GG)", "Bel Ombre (BO)"),
                        levels = c("Gorges (GG)", "Bel Ombre (BO)"))
  )

thresh_df <- data.frame(
  class = factor("Juveniles (0-1 yrs)", levels = levels(surv$class)),
  y     = threshold
)


pop_cols <- c("Gorges (GG)" = "black", "Bel Ombre (BO)" = "darkgreen")
dodge    <- position_dodge(width = 0.6)

surv <- filter(surv, class %in% c("Juveniles (0-1 yrs)", "Pre-breeders (2+ yrs)","Breeders (2+ yrs)")) #subset before plotting (only interested in juv, pre-breeders and breeders)

p_surv <- ggplot(surv, aes(x = Year, y = estimate, colour = Population)) +
  geom_hline(data = thresh_df, aes(yintercept = y), inherit.aes = FALSE,
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4, alpha = 0.3, position = dodge) +
  geom_line(position = dodge) +
  geom_point(position = dodge, size = 1.8) +
  facet_grid(rows = vars(class)) +
  scale_colour_manual(values = pop_cols) +
  scale_x_continuous(breaks = seq(min(surv$Year), max(surv$Year), by = 2)) +
  labs(x = "Year", y = "Survival Estimate (\u03A6)", colour = "Population") +
  theme_bw() +
  theme(
    plot.title   = element_text(color = "#0099f9", size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(color = "black", size = 18, face = "italic"),
    axis.title.y = element_text(color = "black", size = 18, face = "italic"),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y  = element_text(angle = 0, hjust = 1, size = 16),
    strip.text.y = element_text(size = 14),
    legend.position = "top",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )
ggsave(file.path(output_dir, "survival_plot.png"), p_surv,
       width = 13, height = 10, dpi = 300)
print(p_surv)

message("Done. Output files written to '", output_dir, "/'.")
