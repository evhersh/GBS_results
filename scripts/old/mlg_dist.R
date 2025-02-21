# MLG mean dist

# MLG 57 - L41-A_2 looks like it has the lowest mean distance
> mlg57.mat["C27-A_1",]
C27-A_1  C27-A_2  C27-A_3  C27-A_4  C27-A_5  C85-A_1 C85-A_2B  C85-A_3  C85-A_4  C85-A_5  C86-A_1  C86-A_2  C86-A_3 
0.00000 47.57364 47.11736 52.69630 51.45351 59.93318 67.97847 57.69275 56.60709 58.18894 50.67851 51.29771 56.40708 
C86-A_4  C86-A_5  C87-A_1  C87-A_2  C87-A_3  C87-A_4  C87-A_5  C88-A_1  C88-A_2  C88-A_3  C88-A_4  C88-A_5  L41-A_1 
50.85024 49.64753 83.20954 51.52765 49.07059 53.63937 57.07891 51.84097 50.61044 51.31721 51.68633 51.22593 62.20392 
L41-A_2  L41-A_3  L41-A_4  L45-A_2  L45-A_3  L45-A_4  L45-A_5 
50.30240 54.58304 55.70417 56.83994 58.25593 57.31516 56.41229 
> mean(mlg57.mat["C27-A_1",])
[1] 53.362
> mean(mlg57.mat["C27-A_2",])
[1] 52.4475
> mean(mlg57.mat["C27-A_3",])
[1] 53.02056
> mean(mlg57.mat["C27-A_4",])
[1] 53.57804
> mean(mlg57.mat["C27-A_5",])
[1] 53.45413
> mean(mlg57.mat["C85-A_1",])
[1] 59.79365
> mean(mlg57.mat["C85-A_2",])
Error in mlg57.mat["C85-A_2", ] : subscript out of bounds
> mean(mlg57.mat["C85-A_2B",])
[1] 65.58294
> mean(mlg57.mat["C85-A_3",])
[1] 58.22969
> mean(mlg57.mat["C85-A_4",])
[1] 56.83938
> mean(mlg57.mat["C85-A_5",])
[1] 56.30768
> mean(mlg57.mat["C86-A_1",])
[1] 52.83934
> mean(mlg57.mat["C86-A_2",])
[1] 55.36724
> mean(mlg57.mat["C86-A_3",])
[1] 56.26356
> mean(mlg57.mat["C86-A_4",])
[1] 52.36362
> mean(mlg57.mat["C86-A_5",])
[1] 52.90716
> mean(mlg57.mat["C87-A_1",])
[1] 81.54765
> mean(mlg57.mat["C87-A_2",])
[1] 52.61322
> mean(mlg57.mat["C87-A_3",])
[1] 52.33449
> mean(mlg57.mat["C87-A_4",])
[1] 55.54337
> mean(mlg57.mat["C87-A_5",])
[1] 56.59699
> mean(mlg57.mat["C88-A_1",])
[1] 54.03114
> mean(mlg57.mat["C88-A_2",])
[1] 52.88987
> mean(mlg57.mat["C88-A_3",])
[1] 53.35585
> mean(mlg57.mat["C88-A_4",])
[1] 54.18998
> mean(mlg57.mat["C88-A_5",])
[1] 53.15592
> mean(mlg57.mat["L41-A_1",])
[1] 64.1642
> mean(mlg57.mat["L41-A_2",])
[1] 52.06907
> mean(mlg57.mat["L41-A_3",])
[1] 56.30702
> mean(mlg57.mat["L41-A_4",])
[1] 55.80204
> mean(mlg57.mat["L41-A_5",])
Error in mlg57.mat["L41-A_5", ] : subscript out of bounds
> mean(mlg57.mat["L45-A_2",])
[1] 56.35229
> mean(mlg57.mat["L45-A_3",])
[1] 57.1794
> mean(mlg57.mat["L45-A_4",])
[1] 56.28468
> mean(mlg57.mat["L45-A_5",])
[1] 56.07376
> 
  
# MLG 103 - L17-A_3 has lowest mean distance
  
> mlg103.mat["L06-A_1",]
L06-A_1  L06-A_2  L06-A_3  L06-A_4  L06-A_5  L16-A_3  L16-A_5  L17-A_1  L17-A_2  L17-A_3  L17-A_4  L17-A_5 
0.00000 45.16225 46.89868 52.38181 42.10936 70.75739 53.76611 46.10245 48.58835 44.53471 49.47681 45.88489 
> mean(mlg103.mat["L06-A_1",])
[1] 45.4719
> mean(mlg103.mat["L06-A_2",])
[1] 45.79086
> mean(mlg103.mat["L06-A_3",])
[1] 46.24123
> mean(mlg103.mat["L06-A_4",])
[1] 50.63663
> mean(mlg103.mat["L06-A_5",])
[1] 44.76899
> mean(mlg103.mat["L16-A_3",])
[1] 65.03276
> mean(mlg103.mat["L16-A_5",])
[1] 51.13467
> mean(mlg103.mat["L17-A_1",])
[1] 45.70118
> mean(mlg103.mat["L17-A_2",])
[1] 47.75786
> mean(mlg103.mat["L17-A_3",])
[1] 44.32297
> mean(mlg103.mat["L17-A_4",])
[1] 47.63786
> mean(mlg103.mat["L17-A_5",])
[1] 45.65286
> 

# MLG 4 - S03-A_2 has lowest mean dist  
> mlg4.mat["S03-A_1",]
S03-A_1  S03-A_2  S03-A_3  S03-A_4 
0.00000 48.90048 55.16541 51.38150 
> mean(mlg4.mat["S03-A_1",])
[1] 38.86185
> mean(mlg4.mat["S03-A_2",])
[1] 38.22616
> mean(mlg4.mat["S03-A_3",])
[1] 39.94299
> mean(mlg4.mat["S03-A_4",])
[1] 38.29957
> 

# MLG 12 - SM-A_2 has lowest mean dist  
> mlg12.mat["SM-A_1",]
SM-A_1   SM-A_2   SM-A_3   SM-A_4   SM-A_5 
0.00000 63.11237 66.83307 77.59662 64.85332 
> mean(mlg12.mat["SM-A_1",])
[1] 54.47908
> mean(mlg12.mat["SM-A_2",])
[1] 52.94147
> mean(mlg12.mat["SM-A_3",])
[1] 55.4883
> mean(mlg12.mat["SM-A_4",])
[1] 61.91183
> 

# MLG 55 - C43-A_1 has lowest mean distance
> mlg55.mat["C43-A_1",]
C43-A_1  C43-A_2  C43-A_3  C43-A_4  C43-A_5 
0.00000 53.00792 68.15097 81.06642 53.44486 
> mean(mlg55.mat["C43-A_1",])
[1] 51.13403
> mean(mlg55.mat["C43-A_2",])
[1] 51.89453
> mean(mlg55.mat["C43-A_3",])
[1] 56.37829
> mean(mlg55.mat["C43-A_4",])
[1] 64.22175
> mean(mlg55.mat["C43-A_5",])
[1] 52.00055
> 

# MLG 91 - L16-A_2 lowest mean distance  
> mlg91.mat["L16-A_1",]
L16-A_1  L16-A_2  L16-A_4 
0.00000 95.94816 95.11338 
> mean(mlg91.mat["L16-A_1",])
[1] 63.68718
> mean(mlg91.mat["L16-A_2",])
[1] 56.00594

# MLG 104 - L62-A_5 is the only ind from this MLG

# MLG 108 - C23-A_3 has lowest
> mean(mlg108.mat["C23-A_1",])
[1] 47.07106
> mean(mlg108.mat["C23-A_2",])
[1] 43.32903
> mean(mlg108.mat["C23-A_3",])
[1] 41.64993
> mean(mlg108.mat["C23-A_4",])
[1] 43.70445
> mean(mlg108.mat["C23-A_5",])
[1] 43.92722
> 

# MLG 114 - L39-A_2  
> mlg114.mat["L39-A_1",]
L39-A_1  L39-A_2  L39-A_3  L39-A_5 
0.00000 60.17600 68.95904 61.79479 
> mean(mlg114.mat["L39-A_1",])
[1] 47.73246
> mean(mlg114.mat["L39-A_2",])
[1] 44.24506
> mean(mlg114.mat["L39-A_3",])
[1] 49.3135
> mean(mlg114.mat["L39-A_5",])
[1] 45.17815
> 
  
apo.keep.list <- c("L41-A_2", "L17-A_3", "S03-A_2", "SM-A_2", "C43-A_1", "L16-A_2", "L62-A_5", "C23-A_3", "L39-A_2", "L39-A_4")
