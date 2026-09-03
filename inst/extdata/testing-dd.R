# DAS2 examples using 1000 random locations in two dimensions and one
# auxiliary variable


set.seed(2025)

# 1000 points uniformly in [0,1)^2
myLocations <- matrix(runif(1000 * 2), ncol = 2)

# 1000 random directions on the unit circle/sphere
myAuxiliary <- matrix(runif(1000 * 2), ncol = 2)
myAuxiliary <- myAuxiliary / sqrt(rowSums(myAuxiliary^2))


# Draw a spatial balanced (alpha=1) DAS2 sample of n = 20 points
alpha = 1.0;
n = 20;
# candidate samples
sampleMatrix = doubleDAS(myLocations, myAuxiliary, n, alpha);
# first candidate sample
s = sampleMatrix[1,];
s

#######################################

library(spbalDAS)
set.seed(511)
pop <- matrix(runif(400), ncol = 2)
aux <- cbind(pop, runif(nrow(pop)))
dsamp <- doubleDAS(pop, aux, n = 25, alpha = 0.5)
dsamp[1, ]   # one candidate sample of size 25


set.seed(511)
pop <- matrix(runif(4000), ncol = 2)
aux <- matrix(runif(10000), ncol = 5)
samp <- doubleDAS(pop, aux, n = 100, alpha = 0.5)
s = samp[1, ]
s
###########################################

#%% Plot the population (black) and the sample locations (red)
library(ggplot2)

ggplot() +
  geom_point(data = as.data.frame(pop),
             aes(x = V1, y = V2), color = "black", size = 1.5) +
  geom_point(data = as.data.frame(pop[s,]),
             aes(x = V1, y = V2), color = "red", size = 4) +
  labs(title = "Doubly Balanced Sample of n = 100 from 2000 Locations",
       x = "x1 coordinate", y = "x2 coordinate") +
  theme_minimal()








# plotting...
index_matrix <- sampleMatrix
df <- myLocations

library(reshape2)
# Melt the matrix into long format to map indices to colors
index_long <- melt(index_matrix)
colnames(index_long) <- c("row", "col", "index")

# Merge the matrix data with the dataframe coordinates
plot_data <- merge(index_long, df, by.x = "index", by.y = "row.names")

# Create a color column for the row groups of the matrix
plot_data$color_group <- ifelse(plot_data$row == 1, "red", "black")


# Plot the points, colored by the row group of the matrix

ggplot(plot_data, aes(x = V1, y = V2, color = color_group)) +
  scale_color_manual(values=c("black", "red")) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(color = "Matrix Row Group", title = "Spatially Balanced Sample of n = 20 from 1000 Locations",
       subtitle = "100 point pairs, n = 10.",
       caption = "First row colored in RED.") +
  xlab("X Coordinate") +
  ylab("Y Coordinate")


library(ggplot2)

ggplot() +
  geom_point(data = as.data.frame(myLocations),
             aes(x = V1, y = V2), color = "black", size = 1.5) +
  geom_point(data = as.data.frame(myLocations[s, ]),
             aes(x = V1, y = V2), color = "red", size = 4) +
  labs(title = "Spatially Balanced Sample of n = 20 from 1000 Locations",
       x = "x1 coordinate", y = "x2 coordinate") +
  theme_minimal()




