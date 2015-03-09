library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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

p2 <- ggplot(data, aes(y=Distance, x=TreeHeight, fill=algorithm))  +
     stat_summary(fun.data="my_ci_med", size=0.75, aes(colour = algorithm), position = position_dodge(width = 0.3))

p3 <- ggplot(data, aes(y=Distance, x=TreeHeight, fill=algorithm))  +
     stat_summary(fun.data="my_ci_med", size=1.0, position = "jitter") + aes(colour = algorithm)


p4 <- ggplot(data, aes(y=Distance, x=TreeHeight, fill=algorithm)) +
     stat_summary(fun.data="my_ci_med", size = 0.5, position="jitter", aes(group=algorithm, colour = algorithm) ) 


png(filename=paste(args[1],"median.png",sep="."),width = 1980, height = 1200, units = "px")
multiplot(p1, p2, p3, p4, cols=2)
dev.off()
