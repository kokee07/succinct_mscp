# vidca_mscp
Repositorio dedicato al código para Minimum set cover problem trabajado en el proyecto de título de Jorge Delgado.
Existen tres archivos principales, según el formato de trabajo a realizar:

- **VIDCA_EXHAUSTIVE**: Código destinado al procesamiento y comparación entre la solución Exhaustiva, Greedy y Vidca para dataset generado de manera aleatoria según los argumentos.

  	**Ejecución**
./vexhaustive  <n_sets> <MAX_val>
  **n_sets**: el número de subconjuntos a generar.
  **MAX_val**: el máximo valor númerico mayor posible a generar entre 0 y max_val. 
  No es posible modificar el número de elementos por subconjunto actualmente.
  
- **VIDCA_HEURISTIC**: Código destinado al procesamiento y comparación entre Greedy y Vidca para dataset que respeten formato:
    Rep
    #subsets
    sub1
    ..
    subn
    
    **Ejecución**
./vheuristic  <filename>
  **filename**: ruta al archivo con el dataset a utilizar.

    
- **VIDCA_HEURISTIC_2**: Código destinado al procesamiento y comparación entre Greedy y Vidca para dataset que respeten formato:
    #elems #subsets
    cost cardinalidad sub1
    ...
     cost cardinalidad subn
     
     **Ejecución**
     ./vheuristic  <filename>
  **filename**: ruta al archivo con el dataset a utilizar.
