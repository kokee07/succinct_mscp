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

./vheuristic  <file_name>

  **file_name**: ruta al archivo con el dataset a utilizar.


# Flags 

Los archivos tienen 3 flags codeadas para cambiar el comportamiento de salida de los algoritmos:
- **PRINT**: Valor binario 0/1. Permite activar o desactivar mensajes de seguimiento de cada paso del algoritmo.
  	
    **Sugerencia:** Dejar en 1 si es necesario debugear potenciales errores.

- **CHECK**: Valor binario 0/1. Permite ejecutar las rutinas con prefijo **check** que imprimen en detalle la estructura y bitstring utilizada.

    **Sugerencia:** Dejar en 0 durante experimentos, sólo usar para instancias de revisión manual.

- **TEST**: Valor binario 0/1. Permite ejecutar la rutina de test de memoria para corroborar un output correcto de la heurística.

    **Sugerencia:** Dejar en 0 durante experimentos, sólo usar para instancias de revisión manual.