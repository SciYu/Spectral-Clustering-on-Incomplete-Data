# NeurIPS2023-Spectral-Clustering-on-Incomplete-Data

The source code for the NeurIPS'2023 paper titled "**Boosting Spectral Clustering on Incomplete Data via Kernel Correction and Affinity Learning**".

## Introduction

This Github repository contains the implementation of our proposed *Kernel Correction* (KC) algorithm and a series of *self-expressive affinity learning* algorithms (i.e., *KSL-Pp*, *KSL-Sp*, and *AKLSR*), designed to refine kernel estimation and boost spectral clustering, particularly with incomplete data.

- **Kernel Estimation**: Anchored by the kernel matrix's positive semi-definiteness, the KC algorithm provides robust theoretical guarantee for kernel estimation, yielding improved results in kernel-based spectral clustering.

- **Self-expressive Affinity Learning**: The KSL-Pp and KSL-Sp algorithms employ proximal $p$-norm and Schatten $p$-norm penalties within the affinity learning framework to promote sparser affinity matrices. The AKLSR algorithm, an adaptive extension, iteratively learns high-quality kernel and affinity matrices through a joint optimization approach.

## Folders and files

<pre>
./                              - Top directory.
./README.md                     - This readme file.
./main.m                        - Demo of spectral clustering on incomplete data.

|Dataset/                       - Incomplete data
   ./Yale64_0.8miss.m           - Incomplete Yale64 dataset with 80% random missing

|Data_imputation/               - Baselines of data imputation
   ./impute_zero.m              - ZERO Imputation 
   ./impute_mean.m              - MEAN Imputation
   ./impute_knn.m               - kNN Imputation
   ./impute_em.m                - EM Imputation
   ./impute_svt.m               - SVT matrix completion
   ./impute_gr.m                - GR matrix completion
   ./impute_kfmc.m              - KFMC matrix completion

|Distance_calibration/          - Baselines of distance calibration
   ./calibrate_dc.m             - Double-centering algorithm  
   ./calibrate_trf.m            - Triangle fixing algorithm
   ./calibrate_ee.m             - Euclidean embedding algorithm

|Kernel_correction/             - Our proposed kernel correction algorithm 
   ./correct_kernel.m           - Correct kernel matrices (ours)  
   ./correct_distance.m         - Correct Euclidean distance matrices (ours)

|Spectral_clustering/           - Spectral clustering algorithms
   ./cluster_kssc.m             - Kernel sparse subspace clustering
   ./cluster_klsr.m             - Kernel least-squares representation
   ./cluster_ksl_pp.m           - Kernel self-expressive learning with proximal p-norm (ours)
   ./cluster_ksl_sp.m           - Kernel self-expressive learning with Schatten p-norm (ours)
   ./cluster_aklsr.m            - Adaptive kernel least-squares representation (ours)

|Utils/                         - Evaluation files 
   ./distance.m                 - Euclidean distance estimation on incomplete data
   ./generate_x.m               - Perform data imputation
   ./generate_k.m               - Calculate kernel matrices
   ./make_kNN_graph.m           - Calculate a kNN graph
   ./spectral_cluster.m         - Perform spectral clustering algorithms
   ./eval_cluster.m             - Measure clustering performance by [ACC, NMI, PUR, ARI]
   ./eval_error.m               - Measure estimation errors
   ./eval_recall.m              - Measure recall@top-k
   ./statistic_cluster.m        - Statistical results of [ACC, NMI, PUR, ARI, RMSE, RE, Recall]
</pre>


## Citation

If you find this code useful for your research, please use the following BibTeX entry. (to be published)

```
@inproceedings{yu2023cluster,
  title={Boosting Spectral Clustering on Incomplete Data via Kernel Correction and Affinity Learning},
  author={Yu, Fangchen and Zhao, Runze and Shi, Zhan and Lu, Yiwen and Fan, Jicong and Zeng, Yicheng and Mao, Jianfeng and Li, Wenye},
  booktitle={37th Conference on Neural Information Processing Systems},
  year={2023}
}
```

## Contact

If you have any problems or questions, please contact the author: Fangchen Yu (email: fangchenyu@link.cuhk.edu.cn)
