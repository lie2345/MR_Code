# ============================================================================  
# GSMR (Generalized Summary-data-based Mendelian Randomization) Analysis  
# ============================================================================  
# This script performs GSMR analysis to estimate causal effects using summary  
# statistics from GWAS. GSMR extends traditional MR by accounting for linkage  
# disequilibrium between SNPs and implementing pleiotropy outlier detection.  
#  
# Author: [Xuejiao Hou]  
# Date: [2025-03-25]  
# Paper: [ZDHHC4 Influences Obsessive-Compulsive Disorder Risk Through 
# Imaging-Derived Phenotypes: A Mendelian Randomization Study]  
# ============================================================================  

# -----------------------------------------------------------------------------  
# Function to check and install required packages if not already installed  
# -----------------------------------------------------------------------------  
install_if_missing <- function(package) {  
  # Check if package is already installed; if not, attempt to install it  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {  
    # First try installing from CRAN  
    if (!requireNamespace("BiocManager", quietly = TRUE))  
      install.packages("BiocManager")  
    
    tryCatch({  
      install.packages(package)  
    }, error = function(e) {  
      # If CRAN installation fails, try from Bioconductor  
      BiocManager::install(package)  
    })  
    
    library(package, character.only = TRUE)  
    cat(paste0("Package installed and loaded: ", package, "\n"))  
  } else {  
    cat(paste0("Package loaded: ", package, "\n"))  
  }  
}  

# Load required packages for GSMR analysis  
packages <- c("gsmr2", "data.table", "tidyverse", "ggplot2")  
invisible(sapply(packages, install_if_missing))   

# -----------------------------------------------------------------------------  
# Step 1: Load harmonized data from previous MR analysis  
# -----------------------------------------------------------------------------  
# Import the harmonized exposure-outcome dataset created in the main MR analysis  
harm_rt <- fread("harmonise.txt")  

# -----------------------------------------------------------------------------  
# Step 2: Format data for GSMR analysis  
# -----------------------------------------------------------------------------  
# Convert harmonized data from TwoSampleMR format to GSMR format  
# GSMR requires specific columns for the exposure and outcome data  
gsmr_data1 <- data.frame(  
    SNP = harm_rt$SNP,                             # SNP identifier  
    a1 = harm_rt$effect_allele.exposure,           # Effect allele  
    a2 = harm_rt$other_allele.exposure,            # Other allele  
    a1_freq = harm_rt$eaf.exposure,                # Effect allele frequency  
    bzx_n = harm_rt$samplesize.exposure,           # Exposure sample size  
    bzx = harm_rt$beta.exposure,                   # Effect size of SNP on exposure  
    bzx_se = harm_rt$se.exposure,                  # Standard error of SNP-exposure association  
    bzx_pval = harm_rt$pval.exposure,              # P-value of SNP-exposure association  
    bzy = harm_rt$beta.outcome,                    # Effect size of SNP on outcome  
    bzy_n = harm_rt$samplesize.outcome,            # Outcome sample size  
    bzy_se = harm_rt$se.outcome,                   # Standard error of SNP-outcome association  
    bzy_pval = harm_rt$pval.outcome                # P-value of SNP-outcome association  
)  
gsmr_data <- gsmr_data1  

# -----------------------------------------------------------------------------  
# Step 3: Create LD matrix (using identity matrix assuming independence)  
# -----------------------------------------------------------------------------  
# In practice, you should replace this with actual LD data if available  
# Here we assume SNPs are independent (already clumped in the previous analysis)  
# The identity matrix represents no correlation between SNPs  
ldrho <- diag(nrow(gsmr_data))  
colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id <- snpid <- as.character(gsmr_data$SNP)  

# -----------------------------------------------------------------------------  
# Step 4: Extract variables for GSMR analysis  
# -----------------------------------------------------------------------------  
bzx <- gsmr_data$bzx                 # Effect sizes of SNPs on exposure  
bzx_se <- gsmr_data$bzx_se           # Standard errors of the effects  
bzx_pval <- gsmr_data$bzx_pval       # P-values of the SNP-exposure associations  
bzy <- gsmr_data$bzy                 # Effect sizes of SNPs on outcome  
bzy_se <- gsmr_data$bzy_se           # Standard errors of the effects  
bzy_pval <- gsmr_data$bzy_pval       # P-values of the SNP-outcome associations  

# -----------------------------------------------------------------------------  
# Step 5: Set GSMR parameters  
# -----------------------------------------------------------------------------  
# Define the threshold for significant SNPs (typically 5e-8 for GWAS)  
thresh <- 5e-8                         # Genome-wide significance threshold  
n_ref <- 7703                          # Sample size for LD reference panel  
gwas_thresh <- thresh                  # Significance threshold for selecting SNPs  
multi_snps_heidi_thresh <- 0.01        # Significance threshold for HEIDI-outlier test  
nsnps_thresh <- 3                      # Minimum number of SNPs required for analysis  
heidi_outlier_flag <- TRUE             # Flag to enable HEIDI-outlier test (detects pleiotropic outliers)  
ld_r2_thresh <- 0.05                   # LD r² threshold to remove SNPs in high LD  
ld_fdr_thresh <- 0.05                  # FDR threshold for LD pruning  
gsmr2_beta <- 1                        # Parameter for model version (GSMR2)  

# Define exposure and outcome names for plotting  
exposure_name <- "Exposure"            # Replace with actual exposure name  
outcome_name <- "Outcome"              # Replace with actual outcome name  

# -----------------------------------------------------------------------------  
# Step 6: Run GSMR analysis with error handling  
# -----------------------------------------------------------------------------  
# Use tryCatch to handle potential errors during analysis  
tryCatch({  
    # Print debug information  
    print("Debug information before GSMR:")  
    print(paste("Number of SNPs:", length(bzx)))  
    print(paste("Sample size exposure:", harm_rt$samplesize.exposure[1]))  
    print(paste("Sample size outcome:", harm_rt$samplesize.outcome[1]))  
    print(paste("n_ref value:", n_ref))  
    
    # Perform GSMR analysis  
    # This estimates the causal effect while accounting for LD and removing outliers  
    gsmr_results <- gsmr(  
        bzx, bzx_se, bzx_pval,              # Exposure data  
        bzy, bzy_se, bzy_pval,              # Outcome data  
        ldrho, snp_coeff_id, n_ref,         # LD information  
        heidi_outlier_flag,                 # Pleiotropy outlier test flag  
        gwas_thresh,                        # Significance threshold  
        multi_snps_heidi_thresh,            # HEIDI threshold  
        nsnps_thresh,                       # Minimum SNP threshold  
        ld_r2_thresh, ld_fdr_thresh,        # LD pruning parameters  
        gsmr2_beta                          # Model version  
    )  
    
    print("GSMR analysis completed successfully")  

    # -----------------------------------------------------------------------------  
    # Step 7: Calculate effect estimates and confidence intervals  
    # -----------------------------------------------------------------------------  
    # Calculate beta, odds ratio and their 95% confidence intervals  
    beta <- gsmr_results[["bxy"]]                                      # Causal effect estimate  
    beta_lci <- gsmr_results[["bxy"]] - 1.96 * gsmr_results[["bxy_se"]] # Lower CI  
    beta_uci <- gsmr_results[["bxy"]] + 1.96 * gsmr_results[["bxy_se"]] # Upper CI  
    
    # For binary outcomes, convert to odds ratios  
    or <- exp(gsmr_results[["bxy"]])                                      # Odds ratio  
    or_lci <- exp(gsmr_results[["bxy"]] - 1.96 * gsmr_results[["bxy_se"]]) # Lower OR CI  
    or_uci <- exp(gsmr_results[["bxy"]] + 1.96 * gsmr_results[["bxy_se"]]) # Upper OR CI  
    pvalue <- gsmr_results[["bxy_pval"]]                                   # P-value for causal effect  

    # Print results  
    print(paste("Beta (95% CI):", round(beta, 4),   
                "(", round(beta_lci, 4), "-", round(beta_uci, 4), ")"))  
    print(paste("OR (95% CI):", round(or, 4),   
                "(", round(or_lci, 4), "-", round(or_uci, 4), ")"))  
    print(paste("P-value:", format(pvalue, scientific = TRUE, digits = 3)))  

    # Save GSMR results to file  
    capture.output(gsmr_results, file = "gsmr.txt")  

    # -----------------------------------------------------------------------------  
    # Step 8: Create GSMR scatter plot visualization  
    # -----------------------------------------------------------------------------  
    tryCatch({  
        # Get indices of SNPs used in final analysis after filtering  
        # These are the SNPs that passed the HEIDI-outlier test and other filters  
        filtered_index <- gsmr_results$used_index  
        
        # Create PDF file for plot  
        pdf(file = "gsmr.pdf", width = 8, height = 8)  
        
        # Set plot parameters  
        effect_col <- colors()[75]  # Color for effect points  
        
        # Calculate plot limits to ensure all error bars are visible  
        vals <- c(bzx[filtered_index] - bzx_se[filtered_index],   
                  bzx[filtered_index] + bzx_se[filtered_index])  
        xmin <- min(vals)  
        xmax <- max(vals)  
        vals <- c(bzy[filtered_index] - bzy_se[filtered_index],   
                  bzy[filtered_index] + bzy_se[filtered_index])  
        ymin <- min(vals)  
        ymax <- max(vals)  
        
        # Set plot margins  
        par(mar = c(5, 5, 4, 2))  
        
        # Create scatter plot of SNP effects  
        # Each point represents a genetic variant with its effect on exposure and outcome  
        plot(  
            bzx[filtered_index], bzy[filtered_index],   
            pch = 20, cex = 0.8, bty = "n", cex.axis = 1.1, cex.lab = 1.2,  
            col = effect_col, xlim = c(xmin, xmax), ylim = c(ymin, ymax),  
            xlab = bquote(.(exposure_name) ~ (italic(b[zx]))),  
            ylab = bquote(.(outcome_name) ~ (italic(b[zy])))  
        )  
        
        # Add GSMR causal effect line  
        # This line represents the estimated causal effect  
        abline(0, gsmr_results$bxy, lwd = 1.5, lty = 2, col = "dim grey")  
        
        # Add error bars for each SNP  
        nsnps <- length(bzx[filtered_index])  
        for(j in 1:nsnps) {  
            # Horizontal error bars (exposure)  
            xstart <- bzx[filtered_index[j]] - bzx_se[filtered_index[j]]  
            xend <- bzx[filtered_index[j]] + bzx_se[filtered_index[j]]  
            ystart <- bzy[filtered_index[j]]  
            yend <- bzy[filtered_index[j]]  
            segments(xstart, ystart, xend, yend, lwd = 1.5, col = effect_col)  
            
            # Vertical error bars (outcome)  
            xstart <- bzx[filtered_index[j]]  
            xend <- bzx[filtered_index[j]]  
            ystart <- bzy[filtered_index[j]] - bzy_se[filtered_index[j]]  
            yend <- bzy[filtered_index[j]] + bzy_se[filtered_index[j]]  
            segments(xstart, ystart, xend, yend, lwd = 1.5, col = effect_col)  
        }  
        
        # Add title and causal effect information to plot  
        title_text <- paste("GSMR Analysis:", exposure_name, "→", outcome_name)  
        effect_text <- sprintf("Effect (95%% CI): %.3f (%.3f-%.3f), P = %.2e",   
                              beta, beta_lci, beta_uci, pvalue)  
        title(main = title_text, sub = effect_text, cex.main = 1.2, cex.sub = 0.9)  
        
    }, error = function(e) {  
        message("Error in creating GSMR plot: ", e$message)  
    }, finally = {  
        dev.off()  # Close PDF device  
    })  
}, error = function(e) {  
    message("Error in GSMR analysis: ", e$message)  
})  

# -----------------------------------------------------------------------------  
# Step 9: Save formatted results for reporting  
# -----------------------------------------------------------------------------  
# Create a dataframe with formatted results for easy reporting  
results_table <- data.frame(  
    Method = "GSMR",  
    Beta = beta,  
    Beta_LCI = beta_lci,  
    Beta_UCI = beta_uci,  
    OR = or,  
    OR_LCI = or_lci,  
    OR_UCI = or_uci,  
    P_value = pvalue,  
    nSNPs = length(gsmr_results$used_index)  
)  

# Save formatted results  
write.csv(results_table, "gsmr_results.csv", row.names = FALSE)  

# Print completion message  
cat("\nGSMR analysis completed. Results saved to gsmr.txt, gsmr.pdf, and gsmr_results.csv\n")