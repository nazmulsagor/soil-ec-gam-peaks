pkgs <- c("lubridate", "ggplot2", "mgcv", "outliers", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
library(lubridate)
library(ggplot2)
library(mgcv)
library(outliers)
library(patchwork)

main_folder <- NA  # or the data folder path, e.g. "C:/Users/Name/Desktop/clean DATA"

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir)) script_dir <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NULL)
file_arg   <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
file_arg   <- gsub("~+~", " ", file_arg, fixed = TRUE)
if (length(file_arg)) script_dir <- dirname(normalizePath(file_arg))
if (is.null(script_dir)) script_dir <- getwd()
script_dir <- normalizePath(script_dir)

if (is.na(main_folder)) {
  starts <- unique(c(getwd(), script_dir))
  places <- c(file.path(starts, "."), file.path(starts, ".."), file.path(starts, "../.."))
  found  <- places[file.exists(file.path(places, "Treatments.csv"))]
  if (length(found) == 0) stop("Treatments.csv not found. Type the data folder path in main_folder.")
  main_folder <- normalizePath(found[1])
}

site_years <- data.frame(Year     = c(2023,    2024,    2025,     2025,    2026),
                         Location = c("TAREC", "TAREC", "ESAREC", "TAREC", "ESAREC"))

depths <- data.frame(DepthID = 1:9, Depth_cm = seq(10, 90, by = 10))

treatment <- read.csv(file.path(main_folder, "Treatments.csv"), check.names = FALSE)
treatment$`Sensor ID` <- trimws(as.character(treatment$`Sensor ID`))
treatment$Location    <- trimws(treatment$Location)
treatment$Treatment   <- trimws(treatment$Treatment)

pick_location <- "TAREC"
pick_year     <- 2024
pick_sensor   <- "01209A5A"
pick_depth    <- 1
depth_cm      <- pick_depth * 10

read_sensor <- function(file) {
  x <- read.csv(file, fileEncoding = "UTF-8-BOM")
  x$date <- parse_date_time(x$date, orders = c("ymd HMS", "ymd HM", "ymd",
                                                "mdy HMS", "mdy HM", "mdy"))
  x[order(x$date), ]
}

depth_data <- function(x, DepthID) {
  md <- data.frame(date = x$date,
                   SM   = x[, paste0("SM_",   DepthID)],
                   VIC  = x[, paste0("VIC_",  DepthID)],
                   Temp = x[, paste0("Temp_", DepthID, "_avg")])
  na.omit(md)
}

fd_fit <- function(md) {
  VIC.fd  <- diff(md$VIC)
  SM.fd   <- diff(md$SM)
  Temp.fd <- diff(md$Temp)
  glm(VIC.fd ~ SM.fd + Temp.fd - 1, family = gaussian)
}

gam_peak <- function(md, k = 20) {
  b  <- coef(fd_fit(md))
  b0 <- mean(md$VIC) - b[1] * mean(md$SM) - b[2] * mean(md$Temp)
  md$resid <- md$VIC - (b0 + b[1] * md$SM + b[2] * md$Temp)
  md$index <- seq_len(nrow(md))

  gam_m  <- gam(resid ~ s(index, k = k), data = md, method = "REML")
  md$gam <- as.numeric(fitted(gam_m))

  turn <- c(1, which(diff(sign(diff(md$gam))) != 0) + 1, nrow(md))
  x    <- md$gam[turn]
  if (length(x) < 3 || length(x) > 30) return(NULL)

  dx <- dixon.test(x, opposite = (max(x) - mean(x)) < (mean(x) - min(x)))

  list(data = md, peak = turn[which.max(x)],
       q = as.numeric(dx$statistic), p = dx$p.value, edf = summary(gam_m)$edf)
}

legend_inside <- function(x, y, jx = 1, jy = 1) {
  if (packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside", legend.position.inside = c(x, y),
          legend.justification = c(jx, jy))
  } else {
    theme(legend.position = c(x, y), legend.justification = c(jx, jy))
  }
}

results <- data.frame()

for (i in seq_len(nrow(site_years))) {
  folder <- file.path(main_folder, site_years$Year[i], site_years$Location[i])

  for (f in list.files(folder, pattern = "\\.csv$")) {
    x <- read_sensor(file.path(folder, f))

    for (d in 1:9) {
      b <- coef(fd_fit(depth_data(x, d)))
      results <- rbind(results,
                       data.frame(`Sensor ID`       = sub("\\.csv$", "", f),
                                  Depth            = depths$Depth_cm[d],
                                  VWC_Coefficient  = unname(b[1]),
                                  Temp_Coefficient = unname(b[2]),
                                  Year             = site_years$Year[i],
                                  Location         = site_years$Location[i],
                                  check.names      = FALSE))
    }
  }
}

trt <- treatment[, c("Sensor ID", "Year", "Location", "Treatment", "N Treatment Group")]

table_s1 <- merge(results, trt, by = c("Sensor ID", "Year", "Location"), all.x = TRUE)
table_s1 <- table_s1[order(table_s1$Year, table_s1$Location, table_s1$Treatment,
                           table_s1$`Sensor ID`, table_s1$Depth), ]

cat("\nTable S1. FD model coefficients\n")
print(table_s1, row.names = FALSE)

results$SiteYear <- factor(paste(results$Location, results$Year),
                           levels = c("ESAREC 2025", "ESAREC 2026",
                                      "TAREC 2023", "TAREC 2024", "TAREC 2025"))

plot_data <- rbind(data.frame(results[, c("SiteYear", "Location", "Depth")],
                              Coefficient = results$VWC_Coefficient,
                              Parameter   = "VWC coefficient"),
                   data.frame(results[, c("SiteYear", "Location", "Depth")],
                              Coefficient = results$Temp_Coefficient,
                              Parameter   = "Temperature coefficient"))
plot_data$Parameter <- factor(plot_data$Parameter,
                              levels = c("VWC coefficient", "Temperature coefficient"))
plot_data$Depth_f   <- factor(plot_data$Depth)

site_colours <- c("ESAREC 2025" = "#66C2A4", "ESAREC 2026" = "#006D5B",
                  "TAREC 2023"  = "#FDD49E", "TAREC 2024"  = "#F4A259",
                  "TAREC 2025"  = "#B3541E")
location_colours <- c("ESAREC" = "#1B9E77", "TAREC" = "#D95F02")

theme_journal <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.text        = element_text(colour = "black", size = 7),
        axis.title       = element_text(colour = "black", size = 8),
        axis.line        = element_line(colour = "black", linewidth = 0.3),
        axis.ticks       = element_line(colour = "black", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 8, hjust = 0),
        plot.title       = element_text(face = "bold", size = 9),
        panel.spacing    = unit(4, "mm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3.5, "mm"))

set.seed(1)

fig1 <- ggplot(plot_data, aes(x = Depth_f, y = Coefficient, fill = Location)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, linewidth = 0.3,
               outlier.shape = NA, alpha = 0.85) +
  geom_point(aes(colour = Location),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             size = 0.5, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~ Parameter, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = location_colours) +
  scale_colour_manual(values = c("ESAREC" = "#0B5D45", "TAREC" = "#8A3B00")) +
  labs(title = "Fig. 1. FD model coefficients by soil depth, TAREC vs ESAREC",
       x = "Soil depth (cm)", y = "First-difference regression coefficient") +
  theme_journal

print(fig1)

fig2 <- ggplot(plot_data, aes(x = Depth_f, y = Coefficient, fill = SiteYear)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, linewidth = 0.25,
               outlier.size = 0.4, outlier.alpha = 0.5) +
  facet_grid(Parameter ~ Location, scales = "free_y", switch = "y") +
  scale_fill_manual(values = site_colours) +
  labs(title = "Fig. 2. FD model coefficients by soil depth for each site-year",
       x = "Soil depth (cm)", y = NULL) +
  theme_journal +
  theme(strip.placement = "outside",
        strip.text      = element_text(face = "bold", size = 8),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.3)) +
  guides(fill = guide_legend(nrow = 1))

print(fig2)

pvals <- data.frame()

for (i in seq_len(nrow(site_years))) {
  folder <- file.path(main_folder, site_years$Year[i], site_years$Location[i])

  for (f in list.files(folder, pattern = "\\.csv$")) {
    x <- read_sensor(file.path(folder, f))

    for (d in 1:9) {
      p <- summary(fd_fit(depth_data(x, d)))$coefficients[, 4]
      pvals <- rbind(pvals, data.frame(Site = site_years$Location[i],
                                       Year = site_years$Year[i],
                                       VWC  = p[1] < 0.05,
                                       Temp = p[2] < 0.05))
    }
  }
}

vwc  <- tapply(pvals$VWC,  list(pvals$Site, pvals$Year), mean) * 100
temp <- tapply(pvals$Temp, list(pvals$Site, pvals$Year), mean) * 100

table1 <- round(rbind(vwc, temp), 1)
table1 <- data.frame(Site     = rownames(table1),
                     Variable = rep(c("VWC", "Temp"), each = nrow(vwc)),
                     table1, check.names = FALSE, row.names = NULL)
table1 <- table1[order(table1$Site, table1$Variable != "VWC"), ]

cat("\nTable 1. % of FD models with a significant coefficient (p < 0.05)\n")
print(table1, row.names = FALSE)

x  <- read_sensor(file.path(main_folder, pick_year, pick_location, paste0(pick_sensor, ".csv")))
g  <- gam_peak(depth_data(x, pick_depth))
ex <- g$data
pk <- ex[g$peak, ]

site_lab <- paste0(pick_location, " ", pick_year, ", ", depth_cm, " cm depth")

cat(site_lab, "| sensor", pick_sensor, "\n",
    "Peak:", format(pk$date, "%d %b %Y"), "| Dixon Q =", round(g$q, 3),
    "| p =", signif(g$p, 3), "\n")

theme_big <- theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.text         = element_text(colour = "black", size = 10),
        axis.title        = element_text(colour = "black", size = 11),
        axis.line         = element_blank(),
        panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks        = element_line(colour = "black", linewidth = 0.4),
        legend.title      = element_blank(),
        legend.text       = element_text(size = 10),
        legend.background = element_blank(),
        legend.key.width  = unit(7, "mm"),
        plot.title        = element_text(size = 12, face = "bold")) +
  legend_inside(0.99, 0.97)

x_scale <- scale_x_datetime(date_labels = "%b", expand = c(0.01, 0))
no_x    <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
y_room  <- scale_y_continuous(expand = expansion(mult = c(0.05, 0.35)))

p_resid <- ggplot(ex, aes(x = date)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(aes(y = resid, colour = "VIC model residual"), linewidth = 0.35) +
  geom_line(aes(y = gam,   colour = "GAM"), linewidth = 1) +
  geom_point(data = pk, aes(y = gam), shape = 21, fill = "red", colour = "black", size = 2.5) +
  scale_colour_manual(values = c("VIC model residual" = "black", "GAM" = "#1B9E77"),
                      breaks = c("VIC model residual", "GAM")) +
  x_scale + y_room +
  labs(x = paste0("Date (", pick_year, ")"), y = "VIC residuals") +
  theme_big

fig3 <- p_resid + labs(title = paste0("Fig. 3. VIC model residuals and GAM: ", site_lab))
print(fig3)

p_vwc <- ggplot(ex, aes(date, SM)) +
  geom_line(aes(colour = "VWC"), linewidth = 0.4) +
  scale_colour_manual(values = c("VWC" = "#7A0177")) +
  x_scale + y_room + labs(y = "VWC (%)") + theme_big + no_x

p_vic <- ggplot(ex, aes(date, VIC)) +
  geom_line(aes(colour = "VIC"), linewidth = 0.4) +
  scale_colour_manual(values = c("VIC" = "#D95F02")) +
  x_scale + y_room + labs(y = "VIC (-)") + theme_big + no_x

p_temp <- ggplot(ex, aes(date, Temp)) +
  geom_line(aes(colour = "Temperature"), linewidth = 0.4) +
  scale_colour_manual(values = c("Temperature" = "grey40")) +
  x_scale + y_room + labs(y = expression(Temperature~(degree*F))) + theme_big + no_x

fig4 <- p_vwc / p_vic / p_temp / p_resid +
  plot_annotation(title = paste0("Fig. 4. VWC, VIC, temperature and GAM of residuals: ", site_lab),
                  theme = theme(plot.title = element_text(size = 12, face = "bold")))

print(fig4)

peak_table <- data.frame()

for (i in seq_len(nrow(site_years))) {
  folder <- file.path(main_folder, site_years$Year[i], site_years$Location[i])

  for (f in list.files(folder, pattern = "\\.csv$")) {
    x <- read_sensor(file.path(folder, f))

    for (d in 1:9) {
      g <- tryCatch(gam_peak(depth_data(x, d)), error = function(e) NULL)
      if (is.null(g)) next

      pk   <- g$data[g$peak, ]
      type <- "peak"
      if (g$peak == 1)            type <- "series start"
      if (g$peak == nrow(g$data)) type <- "series end"

      peak_table <- rbind(peak_table,
                          data.frame(`Sensor ID`  = sub("\\.csv$", "", f),
                                     Depth       = depths$Depth_cm[d],
                                     k           = 20,
                                     edf         = round(g$edf, 2),
                                     peak.date   = format(pk$date, "%m/%d/%Y %H:%M"),
                                     peak.value  = round(pk$gam, 2),
                                     peak.type   = type,
                                     dixon.q     = round(g$q, 4),
                                     dixon.p     = signif(g$p, 4),
                                     significant = g$p < 0.05 & g$q >= 0.5,
                                     Location    = site_years$Location[i],
                                     Year        = site_years$Year[i],
                                     check.names = FALSE))
    }
  }
}

cat("\nTable 2. Highest GAM peak and Dixon Q test\n")
print(peak_table, row.names = FALSE)

plant <- data.frame(Location = treatment$Location,
                    Year     = treatment$Year,
                    Planting = as.Date(parse_date_time(treatment$`Planting Date`,
                                                       orders = c("mdy", "ymd", "dmy"))))
plant <- aggregate(Planting ~ Location + Year, data = plant, FUN = min)

peaks <- merge(peak_table, plant, by = c("Location", "Year"))
peaks$DAP      <- as.numeric(as.Date(mdy_hm(peaks$peak.date)) - peaks$Planting)
peaks$SiteYear <- paste0(peaks$Location, " (", peaks$Year, ")")

sig     <- subset(peaks, significant)
sig_lab <- "Sig peak (Dixon test: p < 0.05, Q >= 0.5)"
x_max   <- max(160, ceiling(max(peaks$DAP) / 20) * 20)

peak_box <- function(box_fill, box_line) {
  ggplot(peaks, aes(x = DAP, y = Depth, group = Depth)) +
    geom_boxplot(orientation = "y", width = 5, coef = 100,
                 fill = box_fill, colour = box_line, linewidth = 0.4) +
    geom_jitter(width = 0, height = 1.2, colour = "grey35", size = 1.2, alpha = 0.7) +
    geom_point(data = sig, aes(colour = sig_lab), size = 2.5) +
    scale_colour_manual(values = setNames("#D7191C", sig_lab)) +
    scale_y_reverse(breaks = seq(10, 90, 10)) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 20)) +
    facet_wrap(~ SiteYear, ncol = 1, scales = "free_x") +
    labs(title = "Fig. 5. Timing of the highest GAM peak by soil depth",
         x = "Peak timing (days after planting)", y = "Soil depth (cm)") +
    theme_bw(base_size = 12) +
    theme(axis.text          = element_text(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92"),
          strip.background   = element_blank(),
          strip.text         = element_text(face = "bold", size = 13, hjust = 0),
          plot.title         = element_text(face = "bold", size = 14),
          legend.position    = "bottom",
          legend.title       = element_blank())
}

set.seed(1)
print(peak_box("white",   "black"))
set.seed(1)
print(peak_box("#C6DBEF", "#08306B"))
set.seed(1)
print(peak_box("grey85",  "grey20"))

cat("\nPlanting dates used for days after planting\n")
print(plant, row.names = FALSE)

trt <- data.frame(`Sensor ID` = treatment$`Sensor ID`,
                  Year        = treatment$Year,
                  Location    = treatment$Location,
                  Group       = trimws(treatment$`N Treatment Group`),
                  check.names = FALSE)
trt <- subset(trt, paste(Location, Year) %in% paste(site_years$Location, site_years$Year))

peaks <- merge(peak_table, trt, by = c("Sensor ID", "Year", "Location"), all.x = TRUE)

prof <- data.frame()

for (site in c("TAREC", "ESAREC", "All sites")) {

  pk_s  <- if (site == "All sites") peaks else subset(peaks, Location == site)
  trt_s <- if (site == "All sites") trt   else subset(trt,   Location == site)

  for (grp in c("All sensors", sort(unique(trt_s$Group)))) {

    pk_g <- if (grp == "All sensors") pk_s        else subset(pk_s, Group == grp)
    n    <- if (grp == "All sensors") nrow(trt_s) else sum(trt_s$Group == grp)

    s <- aggregate(significant ~ Depth, data = pk_g, FUN = sum)
    prof <- rbind(prof, data.frame(Site = site, Group = grp, Depth = s$Depth,
                                   n_sig = s$significant, n = n))
  }
}

prof$pct   <- 100 * prof$n_sig / prof$n
prof$Lab   <- paste0(prof$Group, " (n = ", prof$n, ")")
prof$lab_x <- ave(prof$pct, prof$Site, prof$Depth, FUN = max) + 3
prof       <- prof[order(prof$Site, prof$Group, prof$Depth), ]

cat("\nFig. 6 data. Sensors with a significant peak by depth\n")
print(prof[, c("Site", "Group", "Depth", "n_sig", "n", "pct")], row.names = FALSE, digits = 3)

x_top <- max(80, ceiling(max(prof$pct) / 10) * 10)

profile_plot <- function(site, letter, cols) {

  d    <- subset(prof, Site == site)
  labs <- unique(d$Lab[order(d$Group != "All sensors")])
  n    <- d$n[d$Group == "All sensors"][1]

  ggplot(d, aes(x = pct, y = Depth, colour = Lab, shape = Lab, linewidth = Lab)) +
    geom_path() +
    geom_point(size = 2.5) +
    geom_text(data = subset(d, Group == "All sensors"),
              aes(x = lab_x, y = Depth, label = paste0(n_sig, "/", n)),
              hjust = 0, size = 3.2, colour = cols[1], inherit.aes = FALSE) +
    scale_colour_manual(values    = setNames(cols[seq_along(labs)], labs), breaks = labs) +
    scale_shape_manual(values     = setNames(c(16, 15, 17)[seq_along(labs)], labs), breaks = labs) +
    scale_linewidth_manual(values = setNames(c(1.2, 0.6, 0.6)[seq_along(labs)], labs), breaks = labs) +
    scale_y_reverse(breaks = seq(10, 90, 10)) +
    scale_x_continuous(limits = c(0, x_top), breaks = seq(0, x_top, 10)) +
    labs(title = paste0("(", letter, ") ", site, " (n = ", n, " sensors)"), x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(axis.text         = element_text(colour = "black"),
          plot.title        = element_text(face = "bold", size = 12),
          legend.title      = element_blank(),
          legend.background = element_blank(),
          legend.key.width  = unit(8, "mm")) +
    legend_inside(0.97, 0.05, 1, 0)
}

p_a <- profile_plot("TAREC",     "a", c("#D7191C", "#7F0000", "#F4A3A3")) + labs(y = "Soil depth (cm)")
p_b <- profile_plot("ESAREC",    "b", c("#2166AC", "#08306B", "#92C5DE")) +
       labs(x = "Sensors with a significant peak (%)")
p_c <- profile_plot("All sites", "c", c("black",   "grey35",  "grey65"))

fig6 <- (p_a | p_b | p_c) +
  plot_annotation(title = "Fig. 6. Depth profile of significant peaks",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

print(fig6)
