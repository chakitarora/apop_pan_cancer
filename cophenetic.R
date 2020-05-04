library("dendextend")

sim2dist <- function(mx) as.dist(1-mx)

J<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Jaccard_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
D<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Dice_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
O1<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Ochiai-1_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
O2<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Ochiai-2_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
RT<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Rogers-Tanimoto_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
RR<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Russell-Rao_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
S<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/Simpson_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram
K<- as.matrix(read.csv("/Users/macbook/Desktop/APOP_PAN_CANCER/jaccard/matrices/kulczynski-2_BPM.csv", sep=",", row.names=1)) %>% sim2dist %>% hclust(method = 'average') %>% as.dendrogram

dend<- dendlist("Dice" = D, "Jaccard" = J, "Ochiai-1" = O1, "Ochiai-2"=O2, "Rogers-Tanimoto"=RT, "Russell-Rao"=RR, "Simpson"=S, "Kulczynski-2"=K)

a=cor.dendlist(dend,method = "baker") #bakers
plot(RT)
library(corrplot)
corrplot(cor.dendlist(dend), "color", "lower")
