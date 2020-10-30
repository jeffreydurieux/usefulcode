# Fri Oct 30 09:51:14 2020
# Author: Jeffrey Durieux, MSc

### soft clustering example using c-means

library(plotly)
library(e1071) #cmeans

data <- iris[,1:2]

# play around with m to change the fuzziness/softness of the clustering
# hard partitioning --> 1.1, soft clustering --> pick a higher number e.g., 2
cm <- cmeans(x = data, centers = 3, m = 1.1)
dff <- data.frame(data, lab = as.factor( cm$cluster), 
                  fuzzy = round(apply(cm$membership, MARGIN = 1, max),3))


#hover over data points to see the fuzziness of the clustering
dff %>% plot_ly(x=~Sepal.Length, y=~Sepal.Width, color = ~lab, text = ~fuzzy)

# check memberships (value between 0 and 1)
cm$membership %>% round(digits = 3)

# 3d plot with fuzzy clustering membership as z axis. 
dff %>% plot_ly(x=~Sepal.Length, y=~Sepal.Width, color = ~lab, z = ~fuzzy)

