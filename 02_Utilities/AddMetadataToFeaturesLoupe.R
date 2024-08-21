### Playground for Loupe Broswer insert meta into Featrues

testCounts <- as.data.frame(MergedSOFCT@assays$RNA$counts)

testMT <- as.data.frame(MergedSOFCT$percent.mt)
colnames(testMT) <- c("MT_Percent")

testMT$MT_Percent <- testMT$MT_Percent * 1000000

testMT <- as.data.frame(t(testMT))


library(Matrix)
library(loupeR)

newCounts <- as.matrix(rbind(testCounts, testMT))

newCountsMTX <- Matrix(newCounts, sparse = T )

loupeR::setup()
create_loupe(count_mat = newCountsMTX,
             clusters = select_clusters(MergedSOFCT),
             projections = select_projections(MergedSOFCT),
             output_name = "MTFeature"
)
