# WeDIV
An Improved K-means Clustering Algorithm Involving new Similarity Metric EP_dis and a Novel Internal Validation Index.
1) MATLAB is the tool of WeDIV
2) Usage

   [Y_optimal,K_optimal,W_optimal] = WeDIV(data,w_step,KList)
   
  Input:
  
    1. Data is the matrix for data set, in which each row is an observation/sample and each column is a variable/feature/attribute.
    
    2. w_step is the the incremental step size for w in EP_dis.
    
    3. KList is the cluster number, which being a constant when the true cluster number was given, otherwise being a vector as 2:NC.
    
    
  Output:
  
    1. Y_optimal is the sample label after clustering.
    
    2. K_optimal is the optimal number of clusters.
    
    3. W_optimal is the optimal of w in EP_dis.
    
