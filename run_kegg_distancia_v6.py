# run_kegg_distancia.py
"""
Estrategia A: Explorar la vía de degradación de PHE y ANT en KEGG
usando un grafo de reacciones del mapa map00624.

VERSIÓN v4 - Estructura para esquema apaisado propio:
- compounds: lista con flags (is_root, is_convergence, is_final, external_link)
- edges: aristas con edge_style (solid/dashed) y edge_direction (forward/lateral)
- external_connections: salidas a otros mapas (Naphthalene, Benzoate, etc.)
- summary: resumen numérico para el reporte
- Lista ordenada de compuestos alcanzados (v3.1)

Requisito: pip install networkx
"""

import requests
import time
import json
import re
from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

try:
    import networkx as nx
except ImportError:
    print("❌ Falta networkx. Instálalo con: pip install networkx")
    raise SystemExit(1)


# ============================================================================
# COFACTORES UBIQUOS A EXCLUIR
# ============================================================================

COFACTORES_UBICUOS = {
    "C00001", "C00080", "C00002", "C00008", "C00020", "C00009", "C00013",
    "C00003", "C00004", "C00005", "C00006", "C00016", "C00015",
    "C00010", "C00255", "C00011", "C00014", "C00007", "C00027",
    "C00025", "C00026", "C00042",
    "C00076", "C00238", "C00244", "C14818",
}


# ============================================================================
# COMPUESTOS PUENTE: conectan map00624 con otros mapas de KEGG
# ============================================================================
# Estos compuestos aparecen en el mapa 00624 y también en otros mapas.
# Sirven para dibujar "cajas externas" en el esquema apaisado.

COMPUESTOS_PUENTE = {
    "C00805": ("map00362", "Benzoate degradation"),      # Salicylate
    "C00090": ("map00362", "Benzoate degradation"),      # Catechol
    "C00230": ("map00362", "Benzoate degradation"),      # Protocatechuate
    "C03012": ("map00626", "Naphthalene degradation"),   # Naphthalene-1,2-diol
    "C01606": ("map00626", "Naphthalene degradation"),   # Phthalate (conecta con 00626)
    "C00251": ("map00362", "Benzoate degradation"),      # Anthranilate (por si aparece)
}


# ============================================================================
# DATACLASSES
# ============================================================================

@dataclass
class EnzymeInfo:
    ec_number: str
    name: str = "—"
    ko: str = "—"
    substrate_id: str = "—"
    substrate_name: str = "—"
    product_id: str = "—"
    product_name: str = "—"
    reaction_id: str = "—"
    position: str = "—"
    distance: int = -1
    is_convergence: bool = False


@dataclass
class PathwayData:
    source: str
    compound: str
    pathway_id: str
    pathway_name: str
    root: Dict = field(default_factory=dict)
    compounds: List[Dict] = field(default_factory=list)
    enzymes: List[Dict] = field(default_factory=list)
    genes: List[Dict] = field(default_factory=list)
    edges: List[Dict] = field(default_factory=list)
    reactions: List[Dict] = field(default_factory=list)
    distances: Dict = field(default_factory=dict)
    final_compounds: List[Dict] = field(default_factory=list)
    convergence_points: List[Dict] = field(default_factory=list)
    external_connections: List[Dict] = field(default_factory=list)
    compounds_reached: List[Dict] = field(default_factory=list)
    summary: Dict = field(default_factory=dict)

    def to_dict(self):
        return {
            "source": self.source,
            "compound": self.compound,
            "pathway_id": self.pathway_id,
            "pathway_name": self.pathway_name,
            "root": self.root,
            "compounds": self.compounds,
            "enzymes": self.enzymes,
            "genes": self.genes,
            "edges": self.edges,
            "reactions": self.reactions,
            "distances": self.distances,
            "final_compounds": self.final_compounds,
            "convergence_points": self.convergence_points,
            "external_connections": self.external_connections,
            "compounds_reached": self.compounds_reached,
            "summary": self.summary,
        }


# ============================================================================
# EXTRACTOR
# ============================================================================

class KEGGDistanciaExtractor:
    BASE_URL = "https://rest.kegg.jp"

    COMPOUNDS = {
        "PHE": "C11422",
        "ANT": "C14315"
    }

    PAH_MAP = "map00624"
    MAX_DISTANCE = 10

    def __init__(self, delay: float = 0.4):
        self.delay = delay
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "KEGG-Distance-Extractor/4.0 (Academic)"
        })
        self._cache_ec_names: Dict[str, str] = {}
        self._cache_ec_kos: Dict[str, List[str]] = {}
        self._cache_cpd_names: Dict[str, str] = {}
        self._reacciones_sustratos: Dict[str, List[str]] = {}
        self._reacciones_productos: Dict[str, List[str]] = {}

    def _rate_limit(self):
        time.sleep(self.delay)

    # ------------------------------------------------------------------
    # 1. REACCIONES DEL MAPA
    # ------------------------------------------------------------------

    def get_reactions_in_map(self) -> List[str]:
        self._rate_limit()
        url = f"{self.BASE_URL}/link/rn/{self.PAH_MAP}"
        reactions = []
        try:
            response = self.session.get(url, timeout=30)
            if response.status_code == 200:
                for line in response.text.strip().split("\n"):
                    if "\t" in line:
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            rn_id = parts[1].replace("rn:", "")
                            if rn_id not in reactions:
                                reactions.append(rn_id)
            print(f"  → {len(reactions)} reacciones encontradas en {self.PAH_MAP}")
        except Exception as e:
            print(f"  ⚠ Error obteniendo reacciones del mapa: {e}")
        return reactions

    # ------------------------------------------------------------------
    # 2. ECUACIÓN DE UNA REACCIÓN
    # ------------------------------------------------------------------

    def get_reaction_equation(self, rn_id: str) -> Tuple[List[str], List[str], Optional[str]]:
        self._rate_limit()
        url = f"{self.BASE_URL}/get/{rn_id}"
        sustratos, productos, ecuacion = [], [], None

        try:
            response = self.session.get(url, timeout=30)
            if response.status_code != 200:
                return sustratos, productos, ecuacion

            for line in response.text.split("\n"):
                if line.startswith("EQUATION"):
                    ecuacion = line.replace("EQUATION", "").strip()
                    if "<=>" in ecuacion:
                        izq, der = ecuacion.split("<=>")
                    elif "=>" in ecuacion:
                        izq, der = ecuacion.split("=>")
                    elif "<=" in ecuacion:
                        izq, der = ecuacion.split("<=")
                    else:
                        continue
                    sustratos = re.findall(r"C\d{5}", izq)
                    productos = re.findall(r"C\d{5}", der)
                    break
        except Exception as e:
            print(f"    ⚠ Error obteniendo ecuación de {rn_id}: {e}")

        return sustratos, productos, ecuacion

    # ------------------------------------------------------------------
    # 3. CONSTRUIR GRAFO
    # ------------------------------------------------------------------

    def construir_grafo(self, reacciones: List[str]) -> nx.DiGraph:
        G = nx.DiGraph()
        print(f"\n[Grafo] Construyendo grafo con {len(reacciones)} reacciones...")
        print(f"  → Excluyendo {len(COFACTORES_UBICUOS)} cofactores ubicuos")

        cofactores_excluidos_count = 0

        for i, rn_id in enumerate(reacciones, 1):
            if i % 20 == 0:
                print(f"  → Procesando reacción {i}/{len(reacciones)}...")

            sustratos, productos, _ = self.get_reaction_equation(rn_id)

            sustratos_limpios = [s for s in sustratos if s not in COFACTORES_UBICUOS]
            productos_limpios = [p for p in productos if p not in COFACTORES_UBICUOS]

            cofactores_excluidos_count += (len(sustratos) - len(sustratos_limpios))
            cofactores_excluidos_count += (len(productos) - len(productos_limpios))

            self._reacciones_sustratos[rn_id] = sustratos_limpios
            self._reacciones_productos[rn_id] = productos_limpios

            for s in sustratos_limpios:
                for p in productos_limpios:
                    if s != p:
                        G.add_edge(s, p, reaction=rn_id)

        print(f"  → Grafo construido: {G.number_of_nodes()} nodos, "
              f"{G.number_of_edges()} aristas")
        print(f"  → Cofactores excluidos: {cofactores_excluidos_count} apariciones")
        return G

    # ------------------------------------------------------------------
    # 4. COMPUESTOS A N PASOS
    # ------------------------------------------------------------------

    def compuestos_a_distancia(self, G: nx.DiGraph,
                                source: str, max_dist: int) -> Dict[int, Set[str]]:
        resultado = {}
        for d in range(1, max_dist + 1):
            try:
                nodos = nx.descendants_at_distance(G, source, d)
                if nodos:
                    resultado[d] = nodos
                    print(f"    distancia {d}: {len(nodos)} compuestos")
                else:
                    break
            except Exception as e:
                print(f"    ⚠ Error en distancia {d}: {e}")
                break
        return resultado

    # ------------------------------------------------------------------
    # 5. POSICIÓN RELATIVA CON LETRAS
    # ------------------------------------------------------------------

    def calcular_posiciones_relativas(self, G: nx.DiGraph, raiz: str,
                                       max_dist: int) -> Dict[str, List[str]]:
        posiciones: Dict[str, List[str]] = defaultdict(list)

        for d in range(1, max_dist + 1):
            try:
                nodos_d = nx.descendants_at_distance(G, raiz, d)
            except Exception:
                break
            if not nodos_d:
                break

            nodos_ordenados = sorted(nodos_d)
            for idx, nodo in enumerate(nodos_ordenados):
                letra = chr(ord('a') + (idx % 26))
                sufijo = f"{d}{letra}"
                if sufijo not in posiciones[nodo]:
                    posiciones[nodo].append(sufijo)

        return dict(posiciones)

    # ------------------------------------------------------------------
    # 6. CLASIFICACIÓN DE ARISTAS
    # ------------------------------------------------------------------

    @staticmethod
    def clasificar_edge_style(ecs: List[str], n_paralelas: int = 1) -> str:
        """
        Determina si la arista debe dibujarse como sólida o punteada.
        - EC con guión (ej: 1.14.13.-) → "dashed" (enzima no completamente caracterizada)
        - Múltiples reacciones entre el mismo par → "dashed"
        - Si no → "solid"
        """
        if not ecs:
            return "dashed"
        for ec in ecs:
            if "-" in ec:
                return "dashed"
        if n_paralelas > 1:
            return "dashed"
        return "solid"

    @staticmethod
    def clasificar_edge_direction(dist_source: int, dist_target: int) -> str:
        """
        - forward: target está a distancia +1 desde source (avance normal)
        - lateral: target está a distancia mucho mayor (salto)
        - reverse: target está a distancia menor (retroceso)
        - same: misma distancia
        """
        if dist_source < 0 or dist_target < 0:
            return "unknown"
        diff = dist_target - dist_source
        if diff == 1:
            return "forward"
        elif diff > 1:
            return "lateral"
        elif diff == 0:
            return "same"
        else:
            return "reverse"

    # ------------------------------------------------------------------
    # 7. LISTA ORDENADA DE COMPUESTOS ALCANZADOS
    # ------------------------------------------------------------------

    def generar_lista_compuestos_alcanzados(self, G: nx.DiGraph,
                                             compound_name: str,
                                             distancias: Dict[int, Set[str]],
                                             compartidos_ids: Set[str]) -> List[Dict]:
        compound_id = self.COMPOUNDS[compound_name]
        posiciones = self.calcular_posiciones_relativas(G, compound_id, self.MAX_DISTANCE)

        lista = []

        # Compuesto padre
        lista.append({
            "code": compound_id,
            "name": self.get_cpd_name(compound_id),
            "distance": 0,
            "letter": "0",
            "is_convergence": compound_id in compartidos_ids,
        })

        for d in sorted(distancias.keys()):
            nodos_d = sorted(distancias[d])
            for idx, cpd in enumerate(nodos_d):
                letra = chr(ord('a') + (idx % 26))
                pos_str = f"{d}{letra}"
                lista.append({
                    "code": cpd,
                    "name": self.get_cpd_name(cpd),
                    "distance": d,
                    "letter": pos_str,
                    "is_convergence": cpd in compartidos_ids,
                })

        lista.sort(key=lambda x: (x["distance"], x["letter"]))
        return lista

    # ------------------------------------------------------------------
    # 8. HELPERS: NOMBRES Y KOs
    # ------------------------------------------------------------------

    def get_ec_name(self, ec: str) -> str:
        if ec in self._cache_ec_names:
            return self._cache_ec_names[ec]

        self._rate_limit()
        url = f"{self.BASE_URL}/get/ec:{ec}"
        nombre = "—"
        try:
            r = self.session.get(url, timeout=30)
            if r.status_code == 200:
                for line in r.text.split("\n"):
                    if line.startswith("NAME"):
                        nombre = line.replace("NAME", "").strip().split(";")[0].strip()
                        break
        except Exception:
            pass

        self._cache_ec_names[ec] = nombre
        return nombre

    def get_ec_kos(self, ec: str) -> List[str]:
        if ec in self._cache_ec_kos:
            return self._cache_ec_kos[ec]

        self._rate_limit()
        kos = []
        try:
            r = self.session.get(f"{self.BASE_URL}/link/ko/ec:{ec}", timeout=30)
            if r.status_code == 200:
                for line in r.text.strip().split("\n"):
                    if "\t" in line:
                        ko = line.split("\t")[1].replace("ko:", "")
                        if ko not in kos:
                            kos.append(ko)
        except Exception:
            pass

        self._cache_ec_kos[ec] = kos
        return kos

    def get_cpd_name(self, cpd_id: str) -> str:
        if cpd_id in self._cache_cpd_names:
            return self._cache_cpd_names[cpd_id]
        if cpd_id == "—" or not cpd_id.startswith("C"):
            return cpd_id

        self._rate_limit()
        url = f"{self.BASE_URL}/get/cpd:{cpd_id}"
        nombre = cpd_id
        try:
            r = self.session.get(url, timeout=30)
            if r.status_code == 200:
                for line in r.text.split("\n"):
                    if line.startswith("NAME"):
                        nombre = line.replace("NAME", "").strip().split(";")[0].strip()
                        break
        except Exception:
            pass

        self._cache_cpd_names[cpd_id] = nombre
        return nombre

    # ------------------------------------------------------------------
    # 9. ENZIMAS DE UNA REACCIÓN
    # ------------------------------------------------------------------

    def get_enzymes_for_reaction(self, rn_id: str) -> List[str]:
        self._rate_limit()
        url = f"{self.BASE_URL}/link/enzyme/rn:{rn_id}"
        enzymes = []
        try:
            response = self.session.get(url, timeout=30)
            if response.status_code == 200:
                for line in response.text.strip().split("\n"):
                    if "\t" in line:
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            ec = parts[1].replace("ec:", "")
                            if ec not in enzymes:
                                enzymes.append(ec)
        except Exception:
            pass
        return enzymes

    # ------------------------------------------------------------------
    # 10. EXTRACCIÓN POR COMPUESTO
    # ------------------------------------------------------------------

    def extract_for_compound(self, compound_name: str,
                              G: nx.DiGraph,
                              distancias: Dict[int, Set[str]],
                              nodos_otro: Set[str],
                              compartidos_ids: Set[str]) -> PathwayData:
        compound_id = self.COMPOUNDS[compound_name]
        print(f"\n[KEGG-Distancia] Procesando {compound_name} ({compound_id})")

        # Recolectar compuestos relevantes
        todos_compuestos = {compound_id}
        for d, nodos in distancias.items():
            todos_compuestos.update(nodos)
        print(f"  → Total de compuestos relevantes: {len(todos_compuestos)}")

        # Posiciones relativas
        posiciones = self.calcular_posiciones_relativas(
            G, compound_id, self.MAX_DISTANCE
        )

        # Mapa de distancias inverso (compuesto -> distancia)
        cpd_a_distancia = {compound_id: 0}
        for d, nodos in distancias.items():
            for n in nodos:
                if n not in cpd_a_distancia:
                    cpd_a_distancia[n] = d

        # Reacciones relevantes
        reacciones_relevantes: Set[str] = set()
        for cpd in todos_compuestos:
            for _, _, data in G.out_edges(cpd, data=True):
                if "reaction" in data:
                    reacciones_relevantes.add(data["reaction"])
            for _, _, data in G.in_edges(cpd, data=True):
                if "reaction" in data:
                    reacciones_relevantes.add(data["reaction"])

        print(f"  → {len(reacciones_relevantes)} reacciones relevantes")

        # Contar reacciones paralelas entre pares (sustrato, producto)
        pares_count = defaultdict(int)
        for rn in reacciones_relevantes:
            sustratos = self._reacciones_sustratos.get(rn, [])
            productos = self._reacciones_productos.get(rn, [])
            for s in sustratos:
                for p in productos:
                    if s != p:
                        pares_count[(s, p)] += 1

        # Construir info enriquecida
        enzimas_info: Dict[str, EnzymeInfo] = {}
        reacciones_info = []
        edges_info: List[Dict] = []

        for i, rn in enumerate(sorted(reacciones_relevantes), 1):
            if i % 20 == 0:
                print(f"    → Procesando reacción {i}/{len(reacciones_relevantes)}...")

            ecs = self.get_enzymes_for_reaction(rn)
            sustratos = self._reacciones_sustratos.get(rn, [])
            productos = self._reacciones_productos.get(rn, [])

            sust_principal = next((s for s in sustratos if s in todos_compuestos), None)
            prod_principal = next((p for p in productos if p in todos_compuestos), None)

            # Info de reacción
            reacciones_info.append({
                "reaction_id": rn,
                "enzymes": ecs,
                "substrates": sustratos,
                "products": productos,
                "main_substrate": sust_principal,
                "main_product": prod_principal
            })

            # Info de aristas (una por cada par sustrato -> producto)
            if sust_principal and prod_principal and sust_principal != prod_principal:
                dist_s = cpd_a_distancia.get(sust_principal, -1)
                dist_p = cpd_a_distancia.get(prod_principal, -1)
                n_paralelas = pares_count.get((sust_principal, prod_principal), 1)

                # Tomar el primer EC como representativo (o "—")
                ec_repr = ecs[0] if ecs else "—"
                ec_name_repr = self.get_ec_name(ec_repr) if ec_repr != "—" else "—"
                ko_repr = ", ".join(self.get_ec_kos(ec_repr)[:3]) if ec_repr != "—" else "—"

                # External link del producto (si aplica)
                external_link = None
                if prod_principal in COMPUESTOS_PUENTE:
                    map_id, map_name = COMPUESTOS_PUENTE[prod_principal]
                    external_link = {"map": map_id, "name": map_name}

                edges_info.append({
                    "from_id": sust_principal,
                    "from_name": self.get_cpd_name(sust_principal),
                    "to_id": prod_principal,
                    "to_name": self.get_cpd_name(prod_principal),
                    "reaction_id": rn,
                    "ec_number": ec_repr,
                    "ec_name": ec_name_repr,
                    "ko": ko_repr,
                    "edge_style": self.clasificar_edge_style(ecs, n_paralelas),
                    "edge_direction": self.clasificar_edge_direction(dist_s, dist_p),
                    "distance_source": dist_s,
                    "distance_target": dist_p,
                    "is_convergence_edge": (sust_principal in nodos_otro or
                                             prod_principal in nodos_otro),
                    "external_link": external_link,
                })

            if not ecs:
                continue

            # Distancia del sustrato principal
            dist_sust = cpd_a_distancia.get(sust_principal, -1) if sust_principal else -1

            # Posición relativa
            if sust_principal and sust_principal in posiciones:
                pos_str = ", ".join(posiciones[sust_principal])
            else:
                pos_str = "—"

            # Es convergencia?
            es_conv = (
                (sust_principal is not None and sust_principal in nodos_otro) or
                (prod_principal is not None and prod_principal in nodos_otro)
            )

            for ec in ecs:
                if ec not in enzimas_info:
                    enzimas_info[ec] = EnzymeInfo(
                        ec_number=ec,
                        name=self.get_ec_name(ec),
                        ko=", ".join(self.get_ec_kos(ec)[:3]) if self.get_ec_kos(ec) else "—",
                        substrate_id=sust_principal or "—",
                        substrate_name=self.get_cpd_name(sust_principal) if sust_principal else "—",
                        product_id=prod_principal or "—",
                        product_name=self.get_cpd_name(prod_principal) if prod_principal else "—",
                        reaction_id=rn,
                        position=pos_str,
                        distance=dist_sust,
                        is_convergence=es_conv
                    )

        print(f"  → {len(enzimas_info)} enzimas únicas")
        print(f"  → {len(edges_info)} aristas con info de estilo")

        # Compuestos finales (nodos hoja)
        finales_ids = []
        for cpd in todos_compuestos:
            if cpd == compound_id:
                continue
            out_edges = list(G.out_edges(cpd))
            if not out_edges:
                finales_ids.append(cpd)

        finales = []
        for cpd in finales_ids:
            ext = COMPUESTOS_PUENTE.get(cpd)
            finales.append({
                "id": cpd,
                "name": self.get_cpd_name(cpd),
                "distance": cpd_a_distancia.get(cpd, -1),
                "external_link": {"map": ext[0], "name": ext[1]} if ext else None,
            })

        # Puntos de convergencia
        convergencia_ids = todos_compuestos & nodos_otro
        convergencia_ids.discard(compound_id)
        convergencia = []
        for cid in sorted(convergencia_ids):
            dist_phe = cpd_a_distancia.get(cid, -1) if compound_name == "PHE" else -1
            dist_ant = cpd_a_distancia.get(cid, -1) if compound_name == "ANT" else -1
            convergencia.append({
                "id": cid,
                "name": self.get_cpd_name(cid),
                "distance_from_this_compound": cpd_a_distancia.get(cid, -1),
            })

        # Conexiones externas
        external_connections = []
        for cpd in todos_compuestos:
            if cpd in COMPUESTOS_PUENTE:
                map_id, map_name = COMPUESTOS_PUENTE[cpd]
                external_connections.append({
                    "compound_id": cpd,
                    "compound_name": self.get_cpd_name(cpd),
                    "target_map": map_id,
                    "target_map_name": map_name,
                    "distance": cpd_a_distancia.get(cpd, -1),
                })

        # Distancias serializables
        dist_serial = {}
        for d, nodos in distancias.items():
            dist_serial[str(d)] = sorted(nodos)

        # Lista de compuestos alcanzados
        lista_compuestos = self.generar_lista_compuestos_alcanzados(
            G, compound_name, distancias, compartidos_ids
        )

        # Lista de compuestos con flags (para el esquema)
        compounds_list = []
        for item in lista_compuestos:
            cid = item["code"]
            compounds_list.append({
                "id": cid,
                "name": item["name"],
                "distance": item["distance"],
                "position": item["letter"],
                "is_root": cid == compound_id,
                "is_convergence": cid in compartidos_ids,
                "is_final": cid in finales_ids,
                "external_link": (lambda ext: {"map": ext[0], "name": ext[1]}
                                   if ext else None)(COMPUESTOS_PUENTE.get(cid)),
            })

        # Root
        root = {
            "id": compound_id,
            "name": self.get_cpd_name(compound_id),
            "distance": 0,
            "position": "0",
        }

        # Summary
        summary = {
            "total_compounds": len(compounds_list),
            "total_enzymes": len(enzimas_info),
            "total_edges": len(edges_info),
            "max_distance_reached": max(distancias.keys()) if distancias else 0,
            "n_convergence_points": len(convergencia),
            "n_final_compounds": len(finales),
            "n_external_links": len(external_connections),
        }

        return PathwayData(
            source="KEGG-Esquematico-v4",
            compound=compound_name,
            pathway_id=self.PAH_MAP,
            pathway_name=f"PAH degradation - {compound_name} (esquemático)",
            root=root,
            compounds=compounds_list,
            enzymes=[e.__dict__ for e in enzimas_info.values()],
            genes=[],
            edges=edges_info,
            reactions=reacciones_info,
            distances=dist_serial,
            final_compounds=finales,
            convergence_points=convergencia,
            external_connections=external_connections,
            compounds_reached=lista_compuestos,
            summary=summary,
        )

    # ------------------------------------------------------------------
    # 11. PIPELINE COMPLETO
    # ------------------------------------------------------------------

    def extract_all(self) -> List[PathwayData]:
        reacciones = self.get_reactions_in_map()
        if not reacciones:
            print("❌ No se pudieron obtener reacciones del mapa")
            return []

        G = self.construir_grafo(reacciones)

        nodos_por_compuesto = {}
        distancias_por_compuesto = {}

        for compound_name, compound_id in self.COMPOUNDS.items():
            print(f"\n[Distancia] Compuestos a distancia desde "
                  f"{compound_name} ({compound_id}):")
            distancias = self.compuestos_a_distancia(G, compound_id, self.MAX_DISTANCE)

            todos = {compound_id}
            for d, nodos in distancias.items():
                todos.update(nodos)

            nodos_por_compuesto[compound_name] = todos
            distancias_por_compuesto[compound_name] = distancias

        # Compuestos compartidos entre PHE y ANT
        compartidos_ids = nodos_por_compuesto["PHE"] & nodos_por_compuesto["ANT"]
        compartidos_ids.discard(self.COMPOUNDS["PHE"])
        compartidos_ids.discard(self.COMPOUNDS["ANT"])

        print(f"\n[Convergencia] Compuestos compartidos PHE ∩ ANT: {len(compartidos_ids)}")

        resultados = []
        for compound_name in self.COMPOUNDS.keys():
            otro = "ANT" if compound_name == "PHE" else "PHE"
            nodos_otro = nodos_por_compuesto[otro]

            data = self.extract_for_compound(
                compound_name,
                G,
                distancias_por_compuesto[compound_name],
                nodos_otro,
                compartidos_ids
            )
            resultados.append(data)

        return resultados


# ============ EJECUCIÓN ============
if __name__ == "__main__":
    print("=" * 60)
    print("EXTRACCIÓN POR DISTANCIA v4 (estructura para esquema apaisado)")
    print("=" * 60)

    extractor = KEGGDistanciaExtractor(delay=0.4)
    datos = extractor.extract_all()

    output = "resultados_kegg_distancia.json"
    with open(output, "w", encoding="utf-8") as f:
        json.dump([d.to_dict() for d in datos], f, indent=2,
                  default=str, ensure_ascii=False)

    print("\n" + "=" * 60)
    print("RESUMEN")
    print("=" * 60)
    for d in datos:
        s = d.summary
        print(f"\n{d.compound}:")
        print(f"  Compuestos:            {s['total_compounds']}")
        print(f"  Enzimas:               {s['total_enzymes']}")
        print(f"  Aristas:               {s['total_edges']}")
        print(f"  Distancia máxima:      {s['max_distance_reached']}")
        print(f"  Compuestos finales:    {s['n_final_compounds']}")
        print(f"  Puntos convergencia:   {s['n_convergence_points']}")
        print(f"  Conexiones externas:   {s['n_external_links']}")

        if d.external_connections:
            print(f"\n  Conexiones externas detectadas:")
            for ec in d.external_connections:
                print(f"    • {ec['compound_name']} (dist {ec['distance']}) "
                      f"→ {ec['target_map_name']}")

        if d.edges:
            n_solid = sum(1 for e in d.edges if e["edge_style"] == "solid")
            n_dashed = sum(1 for e in d.edges if e["edge_style"] == "dashed")
            n_forward = sum(1 for e in d.edges if e["edge_direction"] == "forward")
            n_lateral = sum(1 for e in d.edges if e["edge_direction"] == "lateral")
            print(f"\n  Aristas por estilo:")
            print(f"    Solid:   {n_solid}")
            print(f"    Dashed:  {n_dashed}")
            print(f"  Aristas por dirección:")
            print(f"    Forward: {n_forward}")
            print(f"    Lateral: {n_lateral}")

    print(f"\n✓ Resultados guardados en: {output}")