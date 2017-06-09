#FROM KJENS-1
#Experiment 1: Graph built from healthy control count data
#Experiment 2: Graph built from alzheimer's patient data
#Experiment 3: Graph built from parkinson's patient data
#Experiment 4: Graph built from alzheimer's and healthy control patient count data
#Experiment 5: Graph built from parkinson's and healthy control patient count data
#Experiment 6: Graph built from alzheimer's and parkinson's patient count data
#Experiment 7: Graph from experiment 1, with node names permuted
#Experiment 8: Graph from experiment 2, with node names permuted
#Experiment 9: Graph from experiment 3, with node names permuted
#Experiment 10: Graph from experiment 4, with node names permuted
#Experiment 11: Graph from experiment 5, with node names permuted
#Experiment 12: Graph from experiment 6, with node names permuted

#KJens1-Serum.Only
#Experiment 1: Graph built from healthy control count data
#Experiment 2: Graph built from alzheimer's patient data
#Experiment 3: Graph built from parkinson's patient data
#Experiment 4: Graph built from alzheimer's and healthy control patient count data
#Experiment 5: Graph built from parkinson's and healthy control patient count data
#Experiment 6: Graph built from alzheimer's and parkinson's patient count data
#Experiment 7: Graph from experiment 1, with node names permuted
#Experiment 8: Graph from experiment 2, with node names permuted
#Experiment 9: Graph from experiment 3, with node names permuted
#Experiment 10: Graph from experiment 4, with node names permuted
#Experiment 11: Graph from experiment 5, with node names permuted
#Experiment 12: Graph from experiment 6, with node names permuted


#FROM TPATE DATASET
#Experiment 1: Graph built from healthy control count data
#Experiment 2: Graph built from colon cancer patient data
#Experiment 3: Graph built from prostrate cancer patient data
#Experiment 4: Graph built from colon cancer and healthy control patient count data
#Experiment 5: Graph built from prostrate cancer and healthy control patient count data
#Experiment 6: Graph built from all cancer patient count data
#Experiment 7: Graph from experiment 1, with node names permuted
#Experiment 8: Graph from experiment 2, with node names permuted
#Experiment 9: Graph from experiment 3, with node names permuted
#Experiment 10: Graph from experiment 4, with node names permuted
#Experiment 11: Graph from experiment 5, with node names permuted
#Experiment 12: Graph from experiment 6, with node names permuted

graphSamples = as.data.frame(rbind( c('exrna1.AD','KJENS1-Healthy', 'Alzheimer'),
                     c('exrna2.AD','KJENS1-Alzheimer', 'Alzheimer'),
                     c('exrna3.AD','KJENS1-Parkinson','Alzheimer'),
                     c('exrna4.AD','KJENS1-Alzheimer and Healthy','Alzheimer'),
                     c('exrna5.AD','KJENS1-Parkinson and Healthy','Alzheimer'),
                     c('exrna6.AD','KJENS1-Alzheimer and Parkinson', 'Alzheimer'),
                     c('exrna1.AD','KJENS1-Healthy [serum only]','Alzheimer'),
                     c('exrna4.AD','KJENS1-Alzheimer and Healthy[serum only]','Alzheimer'),
                     c('exrna5.AD','KJENS1-Parkinson and Healthy [serum only]','Alzheimer'),
                     c('exrna6.AD','KJENS1-Alzheimer and Parkinson [serum only]','Alzheimer'),
                     c('exrna1.PD','KJENS1-Healthy','Parkinson'),
                     c('exrna2.PD','KJENS1-Alzheimer','Parkinson'),
                     c('exrna3.PD','KJENS1-Parkinson','Parkinson'),
                     c('exrna4.PD','KJENS1-Alzheimer and Healthy','Parkinson'),
                     c('exrna5.PD','KJENS1-Parkinson and Healthy','Parkinson'),
                     c('exrna6.PD', 'KJENS1-Alzheimer and Parkinson','Parkinson'),
                     c('kras1.Colon','TPATE-Healthy', 'Colon Cancer'),
                     c('kras2.Colon','TPATE-Colon Cancer', 'Colon Cancer'),
                     c('kras3.Colon','TPATE-Prostate Cancer', 'Colon Cancer'),
                     c('kras4.Colon','TPATE-Colon and Healthy', 'Colon Cancer'),
                     c('kras5.Colon','TPATE-Prostate and Healthy', 'Colon Cancer'),
                     c('kras6.Colon','TPATE-Prostate,Colon,Pancreatic', 'Colon Cancer'),
                     c('kras1.Prostate','TPATE-Healthy','Prostate Cancer'),
                     c('kras2.Prostate','TPATE-Colon Cancer','Prostate Cancer'),
                     c('kras3.Prostate','TPATE-Prostate Cancer','Prostate Cancer'),
                     c('kras4.Prostate','TPATE-Colon and Healthy','Prostate Cancer'),
                     c('kras5.Prostate','TPATE-Prostate and Healthy','Prostate Cancer'),
                     c('kras6.Prostate','TPATE-Prostate,Colon,Pancreatic','Prostate Cancer'),
                     c('kras1.Pancreatic','TPATE-Healthy','Pancreatic Cancer'),
                     c('kras2.Pancreatic','TPATE-Colon Cancer','Pancreatic Cancer'),
                     c('kras3.Pancreatic','TPATE-Prostate Cancer','Pancreatic Cancer'),
                     c('kras4.Pancreatic','TPATE-Colon and Healthy','Pancreatic Cancer'),
                     c('kras5.Pancreatic','TPATE-Prostate and Healthy','Pancreatic Cancer'),
                     c('kras6.Pancreatic','TPATE-Prostate,Colon,Pancreatic','Pancreatic Cancer')
                     
))

colnames(graphSamples) <- c('ExpPrName', 'Derived from','Applied to' )
