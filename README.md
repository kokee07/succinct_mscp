
Repository dedicated to the processing and comparison of algorithms for the MSC Problem with different input formats. The main objective is to obtain empirical results of the two minimal heuristics implemented in Jorge Delgado's thesis work [1], surpassing the classical and state-of-the-art Greedy Algorithm.

**Authors:** Jorge Delgado and Hector Ferrada jorge.delgado01@alumnos.uach.cl, hferrada@inf.uach.cl 

**Description:** To achieve our goal, our algorithm called "Succinct-SC" makes use of bit representation techniques for storing subsets, while comparing them using bitwise operations, which accelerates the computation of each iteration. By applying a heuristic that prioritizes minimal elements, i.e., those that are less represented among all the subsets, and creating a data structure that allows iteration over the elements rather than necessarily over the original subsets, we can achieve iterations that are several times faster and more efficient in their decisions.

**Compilation**: To compile all three cases, use the command Make to compile by the Makefile included. This repository uses the library BasicCDS.cpp develop by H. Ferrada.

**Use case:**
- **VIDCA_EXHAUSTIVE:** Code for processing and comparing the Exhaustive, Greedy, and Vidca solutions for datasets generated randomly based on the given arguments.

   **Execution:** ./VIDCA_EXHAUSTIVE <n_sets> <MAX_val>

  - **n_sets:** the number of subsets to generate.

  - **MAX_val:** the maximum numeric value that can be generated between 0 and max_val. It is not currently possible to modify the number of elements per subset.

- **VIDCA_HEURISTIC:** Code for processing and comparing the Greedy and Vidca solutions for datasets that follow this format:

    Rep\
    #subsets\
    sub1
    ..\
    subn\

  **Execution:** ./VIDCA_HEURISTIC <file_name>

  - **file_name:** the path to the file containing the dataset to be used.

- **VIDCA_NATURAL:** Code for processing and comparing the Greedy and Vidca solutions for datasets that follow this format:

    Rep\
    #universe #subsets\
    weight element_count sub1\
    ..\
    weight element_count subn\

    **NOTE:** Weight and element_count are two variables that are not used in this case, but may be relevant in other datasets.

  **Execution:** ./VIDCA_NATURAL <file_name>

- **mscp_parallelism:** Code for processing and comparing the Greedy and **Parallel** Vidca solutions for datasets that follow this format:

    Rep\
    #universe #subsets\
    sub1\
    ..\
    subn\
  
  **Execution:** ./runMSCP <file_name> <#threads>
  
  **NOTE:** #threads should be an integer value to define the numbers of threads to be used in the execution.

  
[1]. J. Delgado, H. Ferrada and C. Navarro. An Efficient Approximate Greedy Algorithm to Solve the Minimum Set Cover Problem.
Submitted to Computational Sciences Journal [2024].

