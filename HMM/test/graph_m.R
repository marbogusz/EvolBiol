library(lattice)

args <- commandArgs(trailingOnly = TRUE)
png(filename = args[2], width = 2560, height = 1440, units = "px")
AAset <- read.table(args[1], header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
with(AAset, xyplot(infer ~ real | type , group=type, auto.key = TRUE, 
                   aspect = "xy",
                   panel = function(x,y,...)
                   {
                     panel.abline(a=0,b=1, col.line="yellow", lwd=3)
                     panel.xyplot(x,y,...)
                   },
                   grid=TRUE,
                   main = args[3],
                   type = c("p", "smooth")))
