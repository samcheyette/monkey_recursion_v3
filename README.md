# monkey_recursion_v3
A hierarchical bayesian model implementation.

There are 4 levels of abstraction in this model: 


Theta_Participant ~ Dirichlet($\alpha$ $\Beta$)
Strategies | Theta_Participant ~ Multinomial(Theta_Participant)
Responses | Strategies ~ Multinomial(Strategies)


at the bottom are participant's responses, which are
multinomial draws from strategies. strategies (e.g. OOMM) give
a dirichlet distribution over responses. 