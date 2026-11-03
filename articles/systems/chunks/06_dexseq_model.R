dxdSubset <- DEXSeqDataSet(
countData = lichIsoCountSubset,
sampleData = lihcMetaFull,
design= ~ sample + exon + condition:exon + gender:exon,
featureID = txInfoSubset$transcript_id,
groupID = txInfoSubset$gene_id
)

dxdSubset <- estimateSizeFactors(dxdSubset)
dxdSubset <- estimateDispersions(dxdSubset)

dxdSubset <- testForDEU(
dxdSubset,
fullModel = ~ sample + exon + condition:exon + gender:exon,
reducedModel = ~ sample + exon + gender:exon
)

dxdSubset <- estimateExonFoldChanges(dxdSubset, fitExpToVar = "condition")
dxdResult <- DEXSeqResults(dxdSubset)
