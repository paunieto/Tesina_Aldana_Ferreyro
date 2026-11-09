
# **Métodos Extendidos**

### **Métodos para el Objetivo 1. 

### **Métodos para el Objetivo 2. Extracción y Análisis Automatizado de Rutas Metabólicas desde KEGG.**

  Para el análisis de las vías de degradación de fenantreno (PHE) y antraceno (ANT), se desarrolló un script en Python (denominado run_kegg_distancia.py) que automatiza la consulta, descarga y estructuración de datos desde la base de datos KEGG (Kyoto Encyclopedia of Genes and Genomes) a través de su API REST (Representational State Transfer). El objetivo principal de esta herramienta fue reconstruir un grafo dirigido de reacciones enzimáticas para determinar qué compuestos son alcanzables desde PHE y ANT, y a qué distancia topológica se encuentran.

#### 2.I. Obtención de reacciones y construcción del grafo
  El proceso de extracción se centró en el mapa de degradación de hidrocarburos aromáticos policíclicos (PAHs) de KEGG, identificado con el código map00624. Inicialmente, el script consultó la totalidad de las reacciones asociadas a este mapa (aproximadamente 119 reacciones). Para cada una de ellas, se descargó su ecuación química y se parseó la información para identificar los sustratos y los productos involucrados.
  Con el fin de evitar la inclusión de reacciones pertenecientes al metabolismo central (que podrían conectar artificialmente la vía con todo el metabolismo celular), se aplicó un filtro de exclusión sobre una lista predefinida de 26 cofactores ubicuos (por ejemplo, agua, ATP, NAD+, NADP+, FAD, CoA, entre otros). A partir de las ecuaciones filtradas, se construyó un grafo dirigido utilizando la librería networkx, donde los nodos representan compuestos químicos (identificados por sus códigos KEGG, e.g., C11422) y las aristas representan las reacciones enzimáticas que los conectan (sustrato → producto).
  
#### 2.II.  Exploración topológica por distancia (BFS)
  Una vez construido el grafo, se implementó un algoritmo de búsqueda en anchura (BFS, Breadth-First Search) para calcular la distancia topológica (número de pasos o reacciones) desde los compuestos raíz hacia el resto de los nodos. Los compuestos raíz fueron el fenantreno (C11422) y el antraceno (C14315). Se exploraron distancias de hasta 10 pasos desde cada raíz para capturar la mayor cantidad de intermediarios posibles dentro del mapa.
Para facilitar la interpretación espacial de la vía, a cada compuesto alcanzable se le asignó una posición relativa compuesta por su distancia desde la raíz y una letra que indica el orden de ramificación dentro de esa distancia (e.g., 1a, 1b, 2a). Esto permitió distinguir entre las diferentes sub-vías que emergen de un mismo nodo.

#### 2.III. Enriquecimiento de la información (enzimas, genes y compuestos)
  Para cada una de las reacciones identificadas como relevantes (es decir, aquellas que conectan compuestos alcanzables desde PHE o ANT), el script consultó la API de KEGG para extraer información detallada:

    - Enzimas (EC numbers): Se recuperaron los códigos EC de las enzimas que catalizan cada reacción.
    - Nombres de enzimas y genes (KOs): Se obtuvo el nombre recomendado para cada EC y se mapearon los identificadores de ortólogos (KOs) asociados, lo que permite vincular la enzima con su información genómica.

    - Sustratos y productos: Para cada enzima, se registró el compuesto de origen (sustrato) y el compuesto resultante (producto), incluyendo sus nombres legibles.

    - Clasificación de aristas: Las reacciones fueron clasificadas según su evidencia biológica. Aquellas con enzimas completamente caracterizadas (EC completo) y que representaban la única vía entre dos compuestos fueron marcadas como reacciones principales (líneas sólidas negras en la visualización). Por otro lado, las reacciones con enzimas hipotéticas (EC con guiones, e.g., 1.14.13.-) o que representaban rutas alternativas/paralelas fueron marcadas como reacciones secundarias o alternativas (líneas grises en la visualización).

#### 2.IV. Detección de convergencia y conexiones externas
  Un objetivo clave del análisis fue identificar los puntos de convergencia entre las vías de PHE y ANT. Para ello, el script comparó los conjuntos de compuestos alcanzables desde ambas raíces e identificó:

    - Compuestos compartidos (puntos de convergencia): Aquellos que tienen enzimas propias en ambas vías, lo que indica que las rutas se unen y continúan juntas. Estos nodos se destacan en color naranja en las visualizaciones.

    - Compuestos alcanzables por ambos: Aquellos que son destinos finales compartidos (nodos hoja) pero que no tienen más reacciones enzimáticas dentro del mapa map00624. Estos nodos se destacan en color violeta.

Adicionalmente, se implementó una búsqueda de conexiones externas, identificando compuestos del mapa map00624 que también actúan como intermediarios en otras rutas metabólicas de KEGG. Este análisis permitió detectar puntos de salida hacia mapas como la degradación de naftaleno (map00626) y la degradación de benzoato (map00362). Estas conexiones se representaron visualmente como cajas con bordes redondeados de color verde esmeralda, conectadas por flechas punteadas del mismo tono, indicando que la vía de degradación de PHE y ANT puede redirigirse hacia esas rutas externas una vez alcanzados dichos intermediarios.

#### 2.V. Almacenamiento y estructuración de los datos
  Todos los datos extraídos, incluyendo los compuestos alcanzables, sus distancias, las enzimas asociadas, los sustratos/productos, y las conexiones externas, fueron almacenados de forma estructurada en un archivo JSON (resultados_kegg_distancia.json). Esta estructura de datos fue diseñada para ser consumida directamente por un script secundario (el visualizador), lo que permite separar la lógica de extracción de datos de la lógica de representación gráfica, facilitando la reproducibilidad y la escalabilidad del análisis.


### **Métodos para el Objetivo 3. Algoritmos de retrosíntesis para la conexión de la vía de degradación de PAHs con el metabolismo central de Pseudomonas putida.**
	Para el diseño de la ruta metabólica se utilizó la plataforma Galaxy-SynBioCAD, la cual integra herramientas de retrosíntesis como RetroRules, RetroPath2.0 y RP2paths. El objetivo fue encontrar una ruta enzimática que conectara los metabolitos nativos de Pseudomonas putida KT2440 con el ftalato.
El flujo de trabajo se inició con la preparación de los insumos. Primero, se intentó descargar el modelo metabólico a escala genómica de P. putida KT2440 (iJN1463) utilizando la herramienta "Pick SBML Model" de Galaxy. Sin embargo, esta herramienta presentó un error de conexión (gzip: stdin: not in gzip format) al intentar descargar el archivo desde la base de datos BiGG, lo que impedía continuar. Para sortear este problema, se optó por una descarga manual del modelo desde el sitio web de BiGG (https://bigg.ucsd.edu/static/models/iJN1463.xml) y su posterior carga en Galaxy mediante la herramienta "Upload Data".

	Con el modelo cargado, se generó el archivo "Sink" (metabolitos nativos del citoplasma) usando la herramienta "Sink from SBML" (compartimento c). Paralelamente, se obtuvieron las reglas de reacción mediante la herramienta "RRules Parser", seleccionando las reglas de tipo retro- y los diámetros de 2 a 16, lo que generó un archivo con 234.384 reglas.
Una vez obtenidos los tres insumos (modelo, Sink y reglas), se ejecutó la herramienta RetroPath2.0 en su versión 2.3.0. Se seleccionó la versión 2.3.0 porque la versión más reciente presentó un error de compatibilidad con la herramienta RP2paths del servidor. Como molécula objetivo (Source InChI) se utilizó el ftalato (InChI=1S/C8H6O4/...), con una longitud máxima de ruta de 3 pasos.

	El siguiente paso en el flujo estándar era utilizar la herramienta RP2paths para separar la red de reacciones en rutas individuales. Sin embargo, esta herramienta falló repetidamente con un error Fatal error: Exit code 1, independientemente de la versión de RetroPath2.0 utilizada. Este fallo se debió a un desajuste interno en las versiones de las herramientas del servidor de Galaxy-SynBioCAD. Para evitar que este error detuviera el análisis, se implementó una solución alternativa: se utilizó la herramienta rp2biosensor, que toma directamente la salida de RetroPath2.0 y genera una red de reacciones en formato JSON y una visualización interactiva en HTML, sin necesidad de pasar por RP2paths. Finalmente, la red se visualizó usando la herramienta "Visualize" de Galaxy, lo que permitió identificar las enzimas (números EC) involucradas en las rutas hacia el ftalato.



