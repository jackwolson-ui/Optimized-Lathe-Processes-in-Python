Problem Overview
-Modeled a cutting tool as a 1D fin with conduction and convection
-Fit 𝑘(𝜙) to a quadratic and h(V) to a power law from provided datasets
-Discretized the governing ODE with finite differences; applied Dirichlet base and Robin tip boundary conditions
-Computed temperature distribution, max temperature, and max gradient; performed a parameter sweep and exported results to Excel

Problem summary
Developed a Python model to support cutting-tool selection by predicting temperature distribution along the tool modeled as a one-dimensional cooling fin. Fit thermal conductivity as a function of steel carbon content and convection coefficient as a function of airflow speed using provided datasets, then solved the resulting boundary value problem via finite-difference discretization with mixed boundary conditions. Generated temperature profiles, extracted maximum temperature and temperature gradient metrics, and automated a parameter sweep across all material/airflow combinations with results exported to an Excel table
