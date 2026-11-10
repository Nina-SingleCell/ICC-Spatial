library(data.table)
library(ggplot2)
library(SeuratObject)
library(ggthemes)

setwd("")
metadata_file <- "Data/ST/SEP/All_hex_SEP_CelltypeConstrained.RDS"
sep_sets <- list(`N-SEP` = paste0("N_SEP", 1:12), `M-SEP` = paste0("M_SEP", 1:5))
obj <- readRDS(metadata_file)
meta <- as.data.table(obj@meta.data)[, .(SampleID, PatientID, Type, SEP)]

for (nm in names(meta)) set(meta, j = nm, value = as.character(meta[[nm]]))

meta <- meta[Type == "Tumor" & SEP %in% unlist(sep_sets)]
meta[, Compartment := fifelse(SEP %in% sep_sets[["N-SEP"]], "N-SEP", "M-SEP")]

counts <- meta[, .(n_hex = .N), by = .(SampleID, PatientID, Compartment, SEP)]
totals <- dcast(counts[, .(n_hex = sum(n_hex)), by = .(SampleID, PatientID, Compartment)],
                SampleID + PatientID ~ Compartment, value.var = "n_hex", fill = 0)
samples <- totals[`N-SEP` > 0 & `M-SEP` > 0]
samples <- samples[PatientID %in% samples[, .N, by = PatientID][N >= 2, PatientID]]

composition <- rbindlist(lapply(names(sep_sets), function(comp) {
    d <- merge(CJ(SampleID = samples$SampleID, SEP = sep_sets[[comp]]),
               counts[Compartment == comp, .(SampleID, SEP, n_hex)],
               by = c("SampleID", "SEP"), all.x = TRUE)
    d[is.na(n_hex), n_hex := 0L]
    d[, proportion := n_hex / sum(n_hex), by = SampleID]
    d[, Compartment := comp]
    merge(d, samples[, .(SampleID, PatientID)], by = "SampleID")
}))

js_distance <- function(p, q) {
    m <- (p + q) / 2
    kl <- function(x) { k <- x > 0; sum(x[k] * log2(x[k] / m[k])) }
    sqrt(max(0, (kl(p) + kl(q)) / 2))
}
parts <- list()
for (comp in names(sep_sets)) {
    wide <- dcast(composition[Compartment == comp], SampleID ~ SEP, value.var = "proportion")
    mat <- as.matrix(wide[, sep_sets[[comp]], with = FALSE])
    rownames(mat) <- wide$SampleID
    for (pat in unique(samples$PatientID)) {
        ids <- sort(samples[PatientID == pat, SampleID])
        ij <- combn(ids, 2)
        parts[[length(parts) + 1L]] <- data.table(
            PatientID = pat, Sample1 = ij[1, ], Sample2 = ij[2, ], Compartment = comp,
            Distance = apply(ij, 2, function(z) js_distance(mat[z[1], ], mat[z[2], ])))
    }
}
pairs <- rbindlist(parts)
patient_long <- pairs[, .(Mean_distance = mean(Distance)), by = .(PatientID, Compartment)]
patient <- dcast(patient_long, PatientID ~ Compartment, value.var = "Mean_distance")
patient[, Delta := `N-SEP` - `M-SEP`]
exact_ok <- !any(patient$Delta == 0) && !anyDuplicated(abs(patient$Delta))
p_value <- if (all(patient$Delta == 0)) 1 else wilcox.test(
    patient$`N-SEP`, patient$`M-SEP`, paired = TRUE, alternative = "two.sided", exact = exact_ok)$p.value
statistics <- data.table(n_patients = nrow(patient), n_regions = nrow(samples),
    n_region_pairs = nrow(pairs) / 2, mean_N = mean(patient$`N-SEP`),
    mean_M = mean(patient$`M-SEP`), mean_difference = mean(patient$Delta),
    n_N_greater = sum(patient$Delta > 0), paired_Wilcoxon_P = p_value)


colors <- c(`N-SEP` = "#237C91", `M-SEP` = "#D46A43")
pairs[, Weight := 1 / .N, by = .(PatientID, Compartment)]


grid_x <- seq(0, 1, length.out = 501)
density_data <- pairs[, {
    w <- Weight / sum(Weight)
    y <- vapply(grid_x, function(z) sum(w * (
        dnorm((z - Distance) / .06) + dnorm((z + Distance) / .06) +
        dnorm((z - (2 - Distance)) / .06))) / .06, numeric(1))
    area <- sum(diff(grid_x) * (head(y, -1) + tail(y, -1)) / 2)
    .(Distance = grid_x, Density = y / area)
}, by = Compartment]


p_density <- ggplot(density_data, aes(Distance, Density,color = Compartment, fill = Compartment)) +
    geom_area(position = "identity", alpha = .15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = colors, limits = names(colors)) +
    scale_fill_manual(values = colors, limits = names(colors)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
        x = "Within-patient Jensen-Shannon distance",
        y = "Density",
        subtitle = sprintf("Paired Wilcoxon P = %.4g", p_value)
    ) +
    theme_few() +
    theme(
        legend.title = element_blank(),
        legend.position = "top"
    )

p_density
