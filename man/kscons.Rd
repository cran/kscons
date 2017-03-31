\name{estimate}
\alias{estimate}
\title{estimate}
\usage{
estimate(a, rho, iter, burn, msa)
}
\arguments{
\item{a}{amino acid sequence}
\item{rho}{secondary structure}
\item{iter}{number of MCMC iterations}
\item{burn}{number of burn-in}
\item{msa}{multiple structure alingment matrix}
}
\description{
Predict the knob-socket model and the corresponding residual contact map based on the amino acid sequence, denoted by a, and the secondary structure, denoted by rho. The users also need to give the iteration and burn-in numbers for the MCMC chain. The users could also input the multiple structure alignment information, by incorporating it into a matrix, where 1) the number of rows and columns of the SYMMETRIC and NON-NEGATIVE matrix should be the length of amino acid sequence; 2) each element give the number of times that each pair of residuals is contacted in the multiple structure alignment; 3) elements are separated by white space. The output is a list, where 1) a is the target primary structure; 2) rho is the target secondary structure; 3) socket is the a 3-by-K matrix, where each column is the position of the triplets; 4) gamma_map is a K-dimension binary vector, it is the MAP estimates of gamma; 5) mpv is a K-diemension vector, it is the estimated probability vector of gamma; 6) $contact_map is a L-by-L binary matrix, it is the contact map corresponding to $gamma_map; 7) $mpm is a L-by-L matrix, it is the estimated probability matrix of contact map.
}
\examples{
a <- "GRIAFTDADDVAILTYVKENARSPSSVTGNALWKAMEKSSLTQHSWQSLKDRYLKHLRG";
rho <- "CCCCCCHHHHHHHHHHHHHCCCTTTTTTTTHHHHHHHHCCTTTCCCHHHHHHHHHHCCC";
iter <- 1000;
burn <- 500;
data(msa_example);
msa <- as.matrix(msa_example);
y <- estimate(a, rho, iter, burn);
}
