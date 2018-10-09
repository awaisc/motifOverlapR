# Using the module_overlap (edges) and module_nodes (nodes), we can now 
# construct the Sankey diagram. A custom colour palette will be defined to 
# reflect the gene modules. 
library(networkD3)
library(dplyr)
library(magrittr)
library(ggalluvial)
library(ggplot2)
# # setwd("/media/awais/Awais/TregDataSets/WGCNA Analysis/")
# save(module_nodes, module_overlap, file = "snakeyModuleInformation.Rdata")
load("snakeyModuleInformation.Rdata")

test<- list()
test$nodes <- dplyr::select(module_nodes, "node")%>%set_colnames("name")
test$links <- dplyr::select(module_overlap, c("PI16modules_num", "StimModules_num", "genes"))%>%
  set_colnames(c("source", "target", "value"))


# module_nodes$node <- module_nodes$node%>%as.character()
# module_nodes$num <- module_nodes$num%>%as.integer()
# module_overlap$PI16modules_num <- module_overlap$PI16modules_num%>%as.integer()
# module_overlap$StimModules_num <- module_overlap$StimModules_num%>%as.integer()
# module_overlap$genes <- module_overlap$genes*100

sankeyNetwork(Links = module_overlap, 
              Nodes = module_nodes, 
              Source = "PI16modules_num", 
              Target = "StimModules_num", 
              Value = "genes", 
              NodeID = "node",
              units = "genes", 
              fontSize = 12,
              width = 909,
              height = 1200,
              fontFamily = "Helvetica",
              LinkGroup = "PI16modules",
              nodeWidth = 30,
              iterations = 10,
              nodePadding = 10)

moduleOverlap <- filter(module_overlap, !(PI16modules == "grey" & StimModules == "grey"))
moduleOverlap$PI16modules <- paste0("PI", moduleOverlap$PI16modules)
moduleOverlap$StimModules <- paste0("Stim", moduleOverlap$StimModules)

ggplot(moduleOverlap ,
       aes(y = genes,
           axis1 = PI16modules  , 
           axis2 = StimModules)) +
  geom_alluvium(aes(fill = isSignificant), 
                width = 1/12) +
  geom_stratum(width = 1/12, fill = moduleOverlapFiltered$PI16modules ,
               color = "grey") +
  scale_x_discrete(limits = c("PI16 modules", "Stim modules"), expand = c(.05, .05)) + 
  ggtitle("Modules preserved between the T stimulation and T PI16 datasets")


