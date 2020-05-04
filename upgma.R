library("dendextend")

sim2dist <- function(mx) as.dist(1-mx)

temp = read.csv("/Users/macbook/Desktop/dice_unfav_yes_1.csv", sep=",", row.names=1)
temp1 <- as.matrix(temp)

temp2<-as.matrix(read.csv("/Users/macbook/Desktop/jaccard_unfav_yes_1.csv", sep=",", row.names=1))
temp3<-as.matrix(read.csv("/Users/macbook/Desktop/ochiai_unfav_yes_1.csv", sep=",", row.names=1))

d=as.dist(sim2dist(temp1))
d=hclust(h,method = 'average') #for UPGMA
d=as.dendrogram(d)
plot(d)

# try these three for the same above task:
#d <- temp1 %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram %>% plot
#d <- temp1 %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
#d <- temp1 %>% sim2dist %>% hclust(method = 'average') %>% plot

d2 <- temp2 %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
d3 <- temp3 %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram

d13<-dendlist(d,d3)
d123<- dendlist("dice" = d, "jaccard" = d2, "ochiai" = d3)
#tanglegram(d13)

#all.equal(d123)

#dist.dendlist(d123)

cor.dendlist(d123) #matrix of cophenetic correlations
library(corrplot)
corrplot(cor.dendlist(d123), "pie", "lower")

#cor_cophenetic(d2, d3) #coph corr bw two dendro

#cor_bakers_gamma(d2,d3) #bakers gamma index

Bk_plot(d, d3, main = "CORRECT Bk plot \n(based on dendrograms)") #bk plot
