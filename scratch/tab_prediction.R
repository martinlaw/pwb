librarian::shelf(data.table, kableExtra, here, devtools)
load_all()
set.seed(16)
scaling <- 100
k10.05 <- addSimonToMeta(nsims=1e5, K=10, theta=0.5, theta0=0.5, theta1=0.6)
k10.03 <- addSimonToMeta(nsims=1e5, K=10, theta=0.3, theta0=0.5, theta1=0.6)
k100.05 <- addSimonToMeta(nsims=1e5, K=100, theta=0.5, theta0=0.5, theta1=0.6)
k100.03 <- addSimonToMeta(nsims=1e5, K=100, theta=0.3, theta0=0.5, theta1=0.6)

k10.05$full.ma
tail(k100.05$full.ma)

########### Table 2 #############
individual.tab <- rbind(k10.05$full.ma, k10.03$full.ma)
scale.names <- c("unw",  "pwb",  "mcse")
individual.tab[, scale.names] <- scaling*individual.tab[, scale.names]
studies.coln <- c(rownames(k10.05$full.ma), rownames(k10.03$full.ma))
individual.tab <- cbind(studies.coln, individual.tab)

# FOR TALKS/HTML (remember to change format="html"):
individual.tab <- kbl(individual.tab,
                   format="latex",
                   digits = c(3,3,3,3,0),
                   align="rrrr",
                   caption="Bias, MCSE (both multiplied by 100) and mean sample size for individual studies, meta-analyses of K=10 studies.",
                   booktabs = T,
                   row.names = FALSE,
                   col.names = c(" ", "Bias", "PWB", "MCSE", "Mean N"),
                   label="individual.bias",
                   escape=FALSE)
# FOR TALKS/HTML (remember to change format="html"):
# smaller.kbl_html <- pack_rows(smaller.kbl, "$\\theta=0.5$", 1, 11, bold = TRUE)
# smaller.kbl_html <- pack_rows(smaller.kbl_html, "$\\theta=0.3$", 12, 22, bold = TRUE)
# smaller.kbl_html <- add_header_above(smaller.kbl_html, c(" " = 1, "Bias" = 2, " " = 1, " " = 1))
# save_kable(smaller.kbl_html, file=here("talks", "bias_individual_studies_100_unw.png"))

individual.tab <- pack_rows(individual.tab, "\\boldsymbol{$\\\\theta=0.5$}", 1, 11, escape = FALSE, bold = FALSE)
individual.tab <- pack_rows(individual.tab, "\\boldsymbol{$\\\\theta=0.3$}", 12, 22, escape = FALSE, bold = FALSE)
#individual.tab <- add_header_above(individual.tab, c(" " = 1, "Bias" = 2, " " = 1, " " = 1))
individual.tab <- kable_styling(individual.tab, latex_options = c("hold_position", "scale_down"))
save_kable(individual.tab, file=here("latex_output", "individ.tex"))






########## TABLE 3 ##############
wide.table.l <- list(k10.05$wide, k10.03$wide, k100.05$wide, k100.03$wide)
wide.table.df <- data.frame(data.table::rbindlist(wide.table.l))

scale.cols <- c("mean.pwb", "bias.ex", "mcse.ex", "bias.inc", "mcse.inc",
                "change", "simon.pw.bias", "pred.change")
wide.table.df[, scale.cols] <- scaling*wide.table.df[, scale.cols]
K <- c(10, 10, 100, 100)
wide.table.df <- cbind(K, wide.table.df)

summary.kbl <- kbl(wide.table.df,
                   format="latex",
                   digits = c(1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                   align="r",
                   caption="Detailed bias and MCSE results for meta-analysis (multiplied by 100), where $W_S$ is the weight of the Simon design in the meta-analysis given that the sum of the weights equals 1, i.e., $\\sum W_i=1$. ``Precision'' is the mean precision of Simon design as a factor of the overall mean precision of remaining studies. We define $\\Delta$ as the difference in bias of the meta-analyses (including minus excluding the study).", #, and define $\\vert \\Delta_{inc} \\vert - \\vert \\Delta_{exc} \\vert$ as the difference in the absolute bias values.
                   booktabs = T,
                   col.names = c("$K$", "\\boldsymbol{$\\theta$}", "$\\overline{PWB_i}$", "Bias", "MCSE", "Bias", "MCSE", "$\\Delta$", "$W_S$", "$\\text{PWB}_{S}$", "$\\text{Predicted }\\Delta$"),
                   label="summary.bias",
                   escape=FALSE)
summary.kbl <- add_header_above(summary.kbl,
                                header=c(" "=2,
                                         "Excluding" = 3,
                                         "Including" = 2,
                                         " "=1,
                                         " "=1,
                                         " "=1,
                                         "$W_S(\\\\text{PWB}_{S} - \\\\text{Bias}_{exc})$"=1),
                                escape = FALSE)
summary.kbl <- add_header_above(summary.kbl, c(" "=2, "Meta-analysis" = 6, "Simon" = 3))
summary.kbl <- kable_styling(summary.kbl, latex_options = c("hold_position", "scale_down"))
save_kable(summary.kbl, file=here("latex_output", "bias_metas_wide.tex"))
