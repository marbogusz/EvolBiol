library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

my_ci_med <- function(x) data.frame(
  y=median(x),
  ymin=quantile(x,0.025),
  ymax=quantile(x,0.975)
)

my_ci <- function(x) data.frame(
  y=mean(x)
)


data <- read.table(args[1], header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

p1 <- ggplot(data, aes(y=Distance, x=TreeHeight, fill=algorithm))  +
     stat_summary(fun.data="my_ci_med", size=0.75, aes(colour = algorithm), position = position_dodge(width = 0.3)) + stat_summary(fun.y = median, geom="line", aes(colour=algorithm))

ggsave(file=paste(args[1],"median.png",sep="."))

