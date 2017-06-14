
# The Hierarchical Bayesian Model is set up this way:


$\alpha \text{ ~} Exponential(\lambda)$
$\beta \text{ ~} Dirichlet(1)$
$\theta^{p} \text{ ~} Dirichlet(\alpha \beta)$
$S | \theta^{p} \text{ ~} Multinomial(\theta^{p})$
$R^{i} | S \text{ ~} Multinomial(S^i)$

$\alpha$ and $\beta$ are population-level hyperparameters: $\beta$ is a vector representing the distribution of strategies across a population (monkeys, kids, tsimane); and $\alpha$ controls how uniform the strategies are for each participant relative to the population. I set $\lambda$ to 0.1 representing a slight preference for non-uniformity in strategies. $\theta^{p}$ is a vector representing a distribution over strategies for each participant. Each participant is then a multinomial over strategies given their $\theta^{p}$. Each participant's responses $R^{i}$ are draws from each of the strategeies $S^{i} \in S$.
a