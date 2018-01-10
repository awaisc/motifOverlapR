PromotersAssoicatedWithEnhancers%>%
  as.data.frame()%>%
  dplyr::select(c("seqnames", "start", "end", "enhancer", "correlation", "FDR"))%>%
  separate(col=`enhancer`, into= c("seqnames2", "start2", "end2"))-> InteractionDataFrameCAGEExpression

SortedTable<-cbind.data.frame(InteractionDataFrameCAGEExpression$seqnames,
                              InteractionDataFrameCAGEExpression$start,
                              InteractionDataFrameCAGEExpression$end,
                              InteractionDataFrameCAGEExpression$seqnames2,
                              InteractionDataFrameCAGEExpression$start2,
                              InteractionDataFrameCAGEExpression$end2,
                              "Interaction",
                              InteractionDataFrameCAGEExpression$correlation*100)

write.table(SortedTable, sep = "\t",
            append=FALSE,
            quote=FALSE,
            col.names = FALSE,
            row.names = FALSE,
            file= "/media/awais/NewDrivewho/UnbiasedTFBSPrediction/Unbiased TFBS Prediction/DataFiles/Interactions/Human/EnhancerPromoterAssoication.bedpe"
            )
EnhancerAndCAGEPRomoterInteraction<-makeGenomicInteractionsFromFile(fn = "/media/awais/NewDrivewho/UnbiasedTFBSPrediction/Unbiased TFBS Prediction/DataFiles/Interactions/Human/EnhancerPromoterAssoication.bedpe", 
                                                                    type = "bedpe")

 mcols(EnhancerAndCAGEPRomoterInteraction)<-cbind.data.frame(mcols(EnhancerAndCAGEPRomoterInteraction), 
                                                             "FDR" = InteractionDataFrameCAGEExpression$FDR)
 
 saveRDS(object = EnhancerAndCAGEPRomoterInteraction, "/media/awais/NewDrivewho/motifOverlapR-master/DataFiles/Interactions/Human/EnhancerPromoterAssoicationRObject")
# 
# InteractionTrack<-InteractionTrack(EnhancerAndCAGEPRomoterInteraction)

CageExpressionGenomicIntearctions <- readRDS(file = "/media/awais/NewDrivewho/motifOverlapR-master/DataFiles/Interactions/Human/EnhancerPromoterAssoicationRObject")

IntearctionTrack<-CageExpressionGenomicIntearctions%>%InteractionTrack()

plotTracks(IntearctionTrack, from= 20000000, to= 25000000, chromosome = "chrX")
