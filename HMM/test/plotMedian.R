library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

my_ci <- function(x) data.frame(
  y=mean(x), 
  ymin=mean(x) - 2 * sd(x), 
  ymax=mean(x) + 2 * sd(x)
)

my_ci_med <- function(x) data.frame(
  y=median(x), 
  ymin=quantile(x,0.025), 
  ymax=quantile(x,0.975)
)

dfr <- read.table(args[1], header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

k <- ggplot(dfr, aes(BirthParam, RFDistance))  + geom_point(color="black",size=1) + facet_grid(algorithm ~ . ) +
     stat_summary(fun.data="my_ci_med", color="red", size=0.3) + ggtitle(args[1]) + coord_cartesian(ylim=c(0, 100))+
 	stat_summary(fun.data="my_ci_med", size=3, colour="darkred", geom="text", show_guide = FALSE, vjust=-1.2, hjust=-0.5,  aes( label=round(..y.., digits=1))) +
 	stat_summary(fun.y=sd, size=3, colour="black", geom="text", show_guide = FALSE, vjust=0, hjust=-0.5,  aes( label=round(..y.., digits=1)))

ggsave(k,file=paste(args[1],"median.png",sep="."))
