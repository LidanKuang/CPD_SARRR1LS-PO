SARRR1LSPO is a novel constrained CPD algorithm with temporal shift-invariance and spatial sparsity and orthonormality constraints (shorted as  sARRR1LS-PO). More specifically, four steps are conducted until convergence for each iteration of the proposed algorithm: 1) use rank-R least-squares fit under spatial phase sparsity constraint to update shared spatial maps after phase de-ambiguity; 2) use orthonormalization constraint to minimize the cross-talk between shared spatial maps; 3) update the aggregating mixing matrix using rank-R least-squares fit; 4) utilize shift-invariant rank-1 least-squares on a series of rank-1 matrices reconstructed by each column of the aggregating mixing matrix to update shared time courses, and subject-specific time delays and intensities. 
It can be applied to complex-valued multi-subject fMRI data analysis.
The main function is SARRR1LSPO.m

Reference:
Li-Dan Kuang, Qiu-Hua Lin, Xiao-Feng Gong, Jianming Zhang, Feng Li, and Vince D. Calhoun, Constrained CPD of Complex-Valued Multi-Subject fMRI Data via Alternating Rank-R and Rank-1 Least Squares, submitted to IEEE Transactions on Neural Systems & Rehabilitation Engineering.
