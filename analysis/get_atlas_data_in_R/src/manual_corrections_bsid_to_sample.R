#Manual corrections to the BSID to RF samplename mappings that 
# were done on 2/11/2017 to correct a few errors that occurred due to a 
# combination of FTP server timeout and a failure of my function 'my.map.BSIDtoSampleName'
# to ensure that each RF metadata file is only read once and that each entry in the
# mapping is matched to only the correct RF metadata file.

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS17D01BLO-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[17] <- as.factor('EXR-KJENS17D01BLO-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_B71_GAGTGG_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[17] <- as.factor('sample_B71_GAGTGG_fastq')

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS109D02CSF-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[18] <- as.factor('EXR-KJENS109D02CSF-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_C92_CACCGG_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[18] <- as.factor('sample_C92_CACCGG_fastq')

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS15D09BLO-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[19] <- as.factor('EXR-KJENS15D09BLO-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_B59_TAATCG_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[19] <- as.factor('sample_B59_TAATCG_fastq')

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS117D10CSF-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[341] <- as.factor('EXR-KJENS117D10CSF-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_C1710_ATGAGC_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[341] <- as.factor('sample_C1710_ATGAGC_fastq')

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS120D07CSF-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[342] <- as.factor('EXR-KJENS120D07CSF-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_C207_CATTTT_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[342] <- as.factor('sample_C207_CATTTT_fastq')

levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`),'EXR-KJENS101D10CSF-BS')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`BS ID`[343] <- as.factor('EXR-KJENS101D10CSF-BS')
levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`) <- c(levels( Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`),'sample_C110_ATTCCT_fastq')
Atlas.all.BStoSample$`KJENS1-aSAH_Project-2016-10-13`$`Sample Name`[343] <- as.factor('sample_C110_ATTCCT_fastq')

