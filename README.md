# Black box optimisation using evolution strategies and tabu search
MATLAB implementation of a few variants of the evolution strategy and tabu search optimisation techniques which are tested on the 5D Rana function.

Detailed information, theoretical discussion and experimental results can be found in "ES_TS_report.pdf".

This was done in the context of the University of Cambridge's Practical Optimisation course (Engineering Tripos, Part IIB), under the supervision of Dr. Geoff Parks.


# Evolution Strategy
We test various recombination and mutation operators (including some personal inventions) and choices of parameters, as well as the Covariance matrix adaptation ES (CMA-ES) technique.
Example of experiments can be found in "Experiment_ES.m" and "Experiment_CMA_ES.m".

<img src="https://user-images.githubusercontent.com/71833961/119832047-be114580-bef5-11eb-835f-26f65aa79a77.png" width="400" height="400"> <img src="https://user-images.githubusercontent.com/71833961/119831883-9ae69600-bef5-11eb-8748-3b6e55a0c4b9.png" width="400" height="400"> 




# Tabu Search

We test various choices of parameters, as well as two variations on the idea of Concentric Tabu Search.
Example of experiments can be found in "Experiment_Tabu.m".

<img src="https://user-images.githubusercontent.com/71833961/119831867-9621e200-bef5-11eb-989e-e85629e04ff8.png" width="400" height="400"> <img src="https://user-images.githubusercontent.com/71833961/119831832-8f936a80-bef5-11eb-9516-035f1bd1d2f0.png" width="400" height="400"> 
