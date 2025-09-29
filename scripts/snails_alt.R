#!/usr/bin/env Rscript
## Script name:
##
## Purpose of script: 
## Data summary and statistics of snail size of Anopolis Lafka Ori.
##
## Author: Savvas Paragkamian
##
## Date Created: 2025-10-01
##
library(tidyverse)
library(readxl)
library(multcompView)
library(factoextra)

snail_data <- read_excel("../data/snails.xlsx",sheet="RAW_data") |>
    mutate(sh_sd=Shell_height/Shell_diameter) |>
    dplyr::select(-c(...9,...10))

snail_data_long <- snail_data |>
    pivot_longer(-c(Individual_ID,Station, ALT),
                 names_to="variables",
                 values_to="value")

breaks <- seq(0, max(snail_data_long$ALT) + 200, by = 200)


labels <- paste0(
  format(head(breaks, -1), scientific = FALSE, big.mark = ""),
  " –",
  format(tail(breaks, -1), scientific = FALSE, big.mark = "")
)


snail_data$alt_bin <- cut(snail_data$ALT,
                               breaks = breaks, labels=labels,
                               right = FALSE)

snail_data_long$alt_bin <- cut(snail_data_long$ALT,
                               breaks = breaks, labels=labels,
                               right = FALSE)

variable_summary <- snail_data_long |>
    group_by(alt_bin,variables) |>
    summarise(mean=mean(value),
              sderr=sd(value,na.rm=T) / sqrt(n()),.groups="keep"
              ) |>
    ungroup() 

write_delim(variable_summary, "../results/variable_summary.tsv",delim="\t")

# --------------------------------------------------------
#################### Statistics ##########################
# --------------------------------------------------------

variables <- unique(snail_data_long$variables)

snails_l_s <- snail_data_long |>
    count(variables,alt_bin)

# ANOVA
aov_results <- snail_data_long |>
    group_by(variables) |>
    do(broom::tidy(aov(value ~ alt_bin,data=.)))

write_delim(stats_results_anova,"../results/stats_results_anova.tsv",delim="\t")

####################### Post hoc test against the control ######################

# Run Tukey HSD per variable
tukey_results <- snail_data_long %>%
  group_by(variables) %>%
  do({
    model <- aov(value ~ alt_bin, data = .)
    tuk <- TukeyHSD(model, "alt_bin")
    broom::tidy(tuk)   # tidy tibble output
  })

write_delim(tukey_results,"../results/tukey_results.tsv",delim="\t")

control_pairwise_sig <- tukey_results |>
    filter(adj.p.value < 0.05) 

write_delim(control_pairwise_sig,"../results/tukey_results_sig.tsv",delim="\t")

letters_per_trait <- snail_data_long %>%
  group_by(variables) %>%
  do({
    model <- aov(value ~ alt_bin, data = .)
    tuk   <- TukeyHSD(model, "alt_bin")
    cld   <- multcompLetters4(model, tuk)
    tibble::tibble(
      alt_bin = names(cld$alt_bin$Letters),
      letters = cld$alt_bin$Letters
    )
  })

#letters_per_trait
trait_means <- snail_data_long %>%
  group_by(variables, alt_bin) %>%
  summarise(max_value = max(value), .groups="drop")

plot_data <- left_join(trait_means, letters_per_trait, 
                       by = c("variables", "alt_bin"))

# --------------------------------------------------------
####################### BOX plot ##########################
# --------------------------------------------------------
vars_a <- unique(snail_data_long$variables)


for (i in seq_along(vars_a)){
    #i=1
    batch_data <- snail_data_long |>
        filter(variables==vars_a[i])

    plot_data_b <- plot_data |>
        filter(variables==vars_a[i])

    print(vars_a[i])
    print("box plot")

    print(i)
    
    fig_bar <- ggplot(batch_data,
                      mapping=aes(x = alt_bin,
                          y = value, fill=alt_bin)) + 
                geom_boxplot(width = 0.6,
                         position = "identity")+
                geom_text(data = plot_data_b, aes(x = alt_bin, y = max_value, label = letters),
                          inherit.aes = FALSE, vjust = 0) +
                geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
                #geom_point(
                #         position = "identity")+
                #position_dodge(width = 0.82)) +
                labs(
                     x = "Altitude",
                     y = vars_a[i]) +
                theme_bw() +
                scale_fill_brewer(palette = "Set3") +
                theme(
                      axis.text.x = element_text(angle = 45,
                                                 hjust = 1,
                                                 size = 11)
                )
    
    ggsave(paste0("../figures/snails_",vars_a[i],"_boxplot.png"),
           plot=fig_bar, 
           height = 20, 
           width = 20,
           dpi = 600, 
           units="cm",
           device="png")
}


# --------------------------------------------------------
########################## PCA ###########################
# --------------------------------------------------------

# Select only numeric morphometric variables
snails_num <- snail_data |> 
    dplyr::select(-c(Individual_ID,ALT,alt_bin,Station))

# Standardize variables (important when measurements are on different scales)
snails_scaled <- scale(snails_num)

# Run PCA
pca_res <- prcomp(snails_scaled, center = TRUE, scale. = TRUE)

# Inspect results
summary(pca_res)        # variance explained
pca_res$rotation        # loadings (contribution of variables)

# as dataframe
pca_scores <- as.data.frame(pca_res$x) %>%
  mutate(alt_bin = snail_data$alt_bin)

fig_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(alt_bin))) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(color = "Altitude bin",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)"))+
         theme_bw() +
         scale_color_brewer(palette = "Set3") +
         theme(
               axis.text.x = element_text(angle = 0,
                                          hjust = 1,
                                          size = 11)
         )

ggsave("../figures/snails_pca.png",
       plot=fig_pca, 
       height = 20, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

# Scree plot (variance explained per PC)
png("../figures/pca_scree.png",
    width = 3000,
    height = 2000,
    res=300,
    units="px")
fviz_eig(pca_res)
dev.off()

# Biplot (specimens + variable loadings)
png("../figures/pca_biplot.png",
    width = 3000,
    height = 2000,
    res=300,
    units="px")
fviz_pca_biplot(pca_res,
                geom.ind = "point",
                habillage = snail_data$alt_bin, # group by species if available
                addEllipses = TRUE)
dev.off()

# Contribution of variables
#fviz_pca_var(pca_res)

