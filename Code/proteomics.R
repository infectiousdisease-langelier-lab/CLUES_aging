# Here I analyze the OLINK proteomic data

# Are main-idea commments

## Are details/steps

# Load packages ----

library(dplyr)
library(tibble)
library(limma)
library(tidyr)
library(ggplot2)
library(ggeffects)

# Read in data =======

cleaned_metadata <- read.csv(
  "./cleaned_CLUES_metadata.csv", row.names = 1, check.names = FALSE
)

p_data <- read.csv("all_IFN_pgml.csv", row.names = 1)


# Merge metadata with proteomic data ====
complete_data <- cleaned_metadata %>%
  inner_join(., p_data, by = c("sampleid" = "ID"))

# Quality control and analysis setup ====

complete_data <- complete_data %>%
  dplyr::select(-c("sledaiscore"))

## Replace NaN with 0
complete_data <- complete_data %>% replace(is.na(.), 0) 

## Isolate proteins, remove no data rows, refactor

# data_df <- complete_data %>%
#   dplyr::filter(across(5:ncol(complete_data), ~ . != "No Data")) %>%
#   dplyr::filter(across(5:ncol(complete_data), ~ . != "> ULOQ")) %>%
#   mutate(across(5:ncol(complete_data), as.numeric)) %>% 
#   mutate(across(c(female, raceeth), as.factor))  

# ## View IFNA2 distribution
# ggplot(data_df, aes(x = IFNA2)) +
#   geom_histogram() +
#   ggtitle("IFNA2 distribution")
# 
# ## View distribution
# ggplot(data_df, aes(x = IFNB1)) +
#   geom_histogram() +
#   ggtitle("IFNB1 distribution")
# 
# ggplot(data_df, aes(x = IFNL1)) +
#   geom_histogram() +
#   ggtitle("IFNL1 distribution")



# Regression analysis

## Set up ggplot2 theme

my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=15, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )


## Analyze

analyze_protein_with_correction <- function(data, proteins) {
  # Initialize list to store results
  results <- list()
  summary_table <- data.frame(Protein = character(0), Slope = numeric(0), 
                              Adjusted_P_Value = numeric(0), P = numeric(0),
                              n = numeric(0))
  # Loop over each protein to compute results
  for (protein in proteins) {
    # Fit linear model
    ## Isolate proteins, remove no data rows, refactor
  
    temp_data <- data %>%
      dplyr::select(age, female, raceeth, !!sym(protein)) %>%
      dplyr::filter(!!sym(protein) != "> ULOQ") %>%
      dplyr::filter(!!sym(protein) != "No Data") %>% 
      mutate(!!sym(protein) := as.numeric(!!sym(protein))) %>%
      mutate(across(c(female, raceeth), as.factor)) 
    
    lt_data <- temp_data %>%
      mutate(!!sym(protein) := log10(!!sym(protein) + 1e-6))
    
    formula <- as.formula(paste(protein, "~", "age + female + raceeth"))
    model <- lm(formula, data = lt_data)
    #print(model)
    
    # Extract coefficient and p-value for age
    coef_summary <- summary(model)$coefficients
    slope <- coef_summary["age", "Estimate"]
    p_value <- coef_summary["age", "Pr(>|t|)"]
    n = nrow(lt_data)
    
    # Store the result temporarily
    results[[protein]] <- list(slope = slope, p_value = p_value, plot = NULL, n = n, lt_data = lt_data)
  }
  
  # Collect all p-values for correction
  p_values <- sapply(results, function(res) res$p_value)
  
  # Correct p-values
  corrected_p_values <- p.adjust(p_values, method = "BH")
  
  for (i in seq_along(proteins)) {
    protein <- proteins[i]
    results[[protein]]$corrected_p_value <- corrected_p_values[i]
    
    
    summary_table <- rbind(
      summary_table,
      data.frame(
        Protein = protein,
        Slope = results[[protein]]$slope,
        Adjusted_P_Value = corrected_p_values[i],
        P = p_values[i],
        n = results[[protein]]$n
      )
    )
    lt <- results[[protein]]$lt_data
    print(dim(lt))
    # Generate scatter plot with regression line
    results[[protein]]$plot <- ggplot(lt, aes_string(x = "age", y = protein)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "#5e3c99") +
      labs(
        title = paste(protein),
        subtitle = paste(
          "Slope =", round(results[[protein]]$slope, 2),
          "| p =", signif(results[[protein]]$p_value, 2),
          "| p_adj =", signif(results[[protein]]$corrected_p_value, 2)
        ),
        x = "Age",
        y = expression(log(concentration ~ (mu*g/ml) + epsilon))
      ) +
      my.theme
  }
  
  return(list(results = results, summary_table = summary_table))
}

proteins <- colnames(complete_data)[5:ncol(complete_data)]

a <- analyze_protein_with_correction(complete_data, proteins)

sum <- a$summary_table

volcano <- ggplot(a$summary_table, aes(x = Slope, y = -log10(Adjusted_P_Value))) +
  geom_point(aes(color = ifelse(Adjusted_P_Value < 0.05, ifelse(Slope > 0, "#e66101", "#5e3c99"), "grey")), size = 3) +  # Corrected color assignment
  scale_color_identity() + 
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  ggtitle("Proteins") + xlab("Slope") +
  ylab(expression(-log[10]("P"[adj]))) + my.theme + theme(axis.text.y = element_text(size = 12), legend.position = "none", legend.box = "vertical",
                                         legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm'), legend.margin=margin(0,0,0,0),
                                         legend.box.margin=margin(-10,-10,-10,-10), legend.text = element_text(size=12)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + xlim(c(-0.02, 0.02))


write.csv(x = a$summary_table, "./summary_table.csv")

ggsave(filename = "./IFNA2.svg", plot = a$results$IFNA2$plot, width = 5, height = 5)

ggsave(filename = "./protein_volcano.svg", plot = volcano, width = 5, height = 5)
