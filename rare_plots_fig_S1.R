# rareplots ----
tiff(filename = "figures/rare_plots.tif",
     width = 20, height = 10, 
     units = "cm", 
     pointsize = 12,
     compression = "lzw",
     res = 300, 
     family = "",
     type = "cairo")

par(mfrow = c(1, 2))
rarecurve(
  t(otu_table(euk_data)),
  step = 50,
  cex = 0.5,
  label = F,
  main = 'Eukaryotes',
  ylab = 'AVS'
)
# abline(v = 5000, lty = 2)

rarecurve(
  t(otu_table(bact_dat)),
  step = 50,
  cex = 0.5,
  label = F,
  main = 'Bacteria',
  ylab = 'AVS'
)

dev.off()

# library(ranacapa)
# ggrare(
#   bact_dat,
#   step = 1000,
#   label = NULL,
#   color = NULL,
#   plot = TRUE,
#   parallel = FALSE
# )
# 
