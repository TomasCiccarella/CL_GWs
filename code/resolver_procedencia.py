#!/usr/bin/env python3
"""Resuelve el conflicto de provenance/ al mergear main en la rama de un modelo.

main y cada rama (monomial, t-model, e-model) agregan entradas al final de provenance/numbers.json
y provenance/claims.yaml, así que al hacer `git merge main` en una rama git suele marcar conflicto
en esos dos archivos. La resolución correcta es la unión: lo de la rama más lo nuevo de main.

Uso, desde la raíz del repositorio y con el merge en conflicto:
    git checkout t-model
    git merge main                          # CONFLICT en provenance/...
    python3 code/resolver_procedencia.py
    git add provenance && git commit --no-edit

numbers.json: conserva las claves de la rama y agrega las que sólo están en main; si una clave está
en los dos lados con valores distintos, se detiene (eso no es un conflicto de formato y hay que
mirarlo a mano). claims.yaml: conserva el texto de la rama y agrega al final las afirmaciones de
main cuyo id no está en la rama; en la sección figures: hace lo mismo con las figuras, por file. Al final comprueba que no haya ids repetidos y que todo número
citado por una afirmación exista.
"""

import json
import subprocess
import sys

import yaml


def show(etapa, ruta):
    """Contenido de un archivo en una etapa del índice durante un merge (:2 rama, :3 main). Si git
    ya mezcló ese archivo sin conflicto, las dos versiones son la del árbol de trabajo."""
    r = subprocess.run(["git", "show", f"{etapa}:{ruta}"], capture_output=True, text=True)
    if r.returncode == 0:
        return r.stdout
    if subprocess.run(["git", "rev-parse", "-q", "--verify", "MERGE_HEAD"], capture_output=True).returncode:
        sys.exit("No hay un merge en curso (¿corriste git merge main?)")
    return open(ruta, encoding="utf-8").read()


def separar(txt):
    """claims.yaml -> (texto de 'claims:', texto de los ítems de 'figures:' o '')."""
    partes = txt.split("\nfigures:\n", 1)
    return partes[0], partes[1] if len(partes) > 1 else ""


def unir_bloques(txt_rama, txt_otra, marca):
    """Los ítems de txt_rama más los de txt_otra cuya clave (id o file) no está en la rama.
    Devuelve el texto completo de la unión si marca es la de figuras, y sólo los nuevos si es la
    de afirmaciones (que se agregan al final del texto de la rama)."""
    def items(txt):
        return [marca + b.rstrip("\n") + "\n" for b in ("\n" + txt).split("\n" + marca)[1:]]
    def clave(item):
        return item[len(marca):].split("\n")[0].strip()
    de_rama = items(txt_rama)
    vistos = {clave(i) for i in de_rama}
    nuevos = [i for i in items(txt_otra) if clave(i) not in vistos]
    return "".join(de_rama + nuevos) if "file" in marca else "".join(nuevos)


def main():
    rama = json.loads(show(":2", "provenance/numbers.json"))
    otra = json.loads(show(":3", "provenance/numbers.json"))
    for k, v in otra.items():
        if k not in rama:
            rama[k] = v
        elif rama[k] != v:
            sys.exit(f"numbers.json[{k}] tiene valores distintos en la rama y en main: resolvelo a mano.")
    with open("provenance/numbers.json", "w", encoding="utf-8") as fh:
        json.dump(rama, fh, indent=2, ensure_ascii=False)

    cl_rama, fig_rama = separar(show(":2", "provenance/claims.yaml"))
    cl_otra, fig_otra = separar(show(":3", "provenance/claims.yaml"))
    nuevos = unir_bloques(cl_rama, cl_otra, "  - id: ")
    figs = unir_bloques(fig_rama, fig_otra, "  - file: ")
    txt = cl_rama.rstrip("\n") + "\n" + "".join(nuevos)
    if figs.strip():
        txt += "\nfigures:\n" + figs.lstrip("\n")
    with open("provenance/claims.yaml", "w", encoding="utf-8") as fh:
        fh.write(txt)

    ids = [c["id"] for c in yaml.safe_load(txt)["claims"]]
    faltan = [k for c in yaml.safe_load(txt)["claims"] for k in c.get("numbers", []) if k not in rama]
    if len(ids) != len(set(ids)) or faltan:
        sys.exit(f"Resultado inconsistente: ids repetidos o números faltantes {faltan}")
    print(f"provenance/ resuelto: {len(ids)} afirmaciones ({nuevos.count('  - id: ')} nuevas de main), "
          f"{len(rama)} números. Ahora: git add provenance && git commit --no-edit")


if __name__ == "__main__":
    main()
