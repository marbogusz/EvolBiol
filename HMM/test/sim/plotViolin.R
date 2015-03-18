library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

my_ci_med <- function(x) data.frame(
  y=median(x),
  ymin=quantile(x,0.025),
  ymax=quantile(x,0.975)
)

my_ci <- function(x) data.frame(
  y=mean(x),
  ymin=mean(x) - sd(x),
  ymax=mean(x) + sd(x)
)


data <- read.table(args[1], header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

p1 <- ggplot(data, aes(y=Distance, x=TreeHeight ))  + 
      geom_violin(scale = "width", size = 0.1, aes(fill = algorithm, group=interaction(TreeHeight,algorithm)), position = position_dodge(width = 0.3))+
      stat_summary(fun.y = median, geom="line", size = 0.1,  aes(colour=algorithm),position = position_dodge(width = 0.3))+
      ggtitle(args[1]) + ylab("RF distance to true")


ggsave(file=paste(args[1],"violin.png",sep="."))

